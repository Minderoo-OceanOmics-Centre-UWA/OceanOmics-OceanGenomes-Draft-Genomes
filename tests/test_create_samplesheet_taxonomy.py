"""master_species is not the only source of a sample's taxonomy.

A miss there used to abort the whole run: create_samplesheet.py refuses to emit
taxon_id='unknown' (FCS-GX is handed `--tax-id unknown` and fails only after
MEGAHIT has burned the SUs), and the only remedy offered was to bulk-load the
taxa by hand. The invertebrate runs draw from most of Metazoa, so that does not
scale -- NOVA_260909_LA lost 67 of its samples to it in one go.

These tests pin the NCBI taxdump fallback that fills the gap, the two resolver
flags this pipeline needs that the mitogenome one must not have, the rule that a
name absent from BOTH sources still fails loudly rather than guessing, and the
rule that a congener's row never supplies the taxid when NCBI knows the sample's
own name.
"""

import contextlib
import csv
import importlib.util
import io
import json
import os
import re
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"

sys.path.insert(0, str(BIN))
from taxdump_lineage import TaxdumpLineage  # noqa: E402,F401


def load_create_samplesheet():
    """Import create_samplesheet.py without requiring psycopg2 to be installed."""
    spec = importlib.util.spec_from_file_location(
        "create_samplesheet_taxonomy_under_test", BIN / "create_samplesheet.py"
    )
    module = importlib.util.module_from_spec(spec)
    patches = {}
    try:
        import psycopg2  # noqa: F401
    except ImportError:
        patches["psycopg2"] = types.ModuleType("psycopg2")
    sys.path.insert(0, str(BIN))
    try:
        with mock.patch.dict(sys.modules, patches):
            spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(BIN))
    return module


cs = load_create_samplesheet()


# A miniature dump rooted at Metazoa, covering the shapes the real run hits:
#   Arthropoda > Malacostraca > Decapoda > Caridea (infraorder) > Chaceon (genus)
#     > Chaceon bicolor (species)
#   Echinodermata > Asteroidea (a class used as the whole nominal name)
#   Placozoa (a phylum with no class rank) > Trichoplax
#   Streptophyta > Morus, outside Metazoa, to prove the clade bound holds
NODES = [
    (1, 1, "no rank"),
    (33208, 1, "kingdom"),
    (100, 33208, "phylum"),      # Arthropoda
    (101, 100, "class"),         # Malacostraca
    (102, 101, "order"),         # Decapoda
    (103, 102, "infraorder"),    # Caridea
    (104, 103, "genus"),         # Chaceon
    (105, 104, "species"),       # Chaceon bicolor
    (110, 33208, "phylum"),      # Echinodermata
    (111, 110, "class"),         # Asteroidea
    (200, 33208, "phylum"),      # Placozoa
    (201, 200, "genus"),         # Trichoplax
    (300, 1, "phylum"),          # Streptophyta
    (301, 300, "genus"),         # Morus, the mulberry
]

NAMES = [
    (1, "root", "scientific name"),
    (33208, "Metazoa", "scientific name"),
    (100, "Arthropoda", "scientific name"),
    (101, "Malacostraca", "scientific name"),
    (102, "Decapoda", "scientific name"),
    (103, "Caridea", "scientific name"),
    (104, "Chaceon", "scientific name"),
    (105, "Chaceon bicolor", "scientific name"),
    (110, "Echinodermata", "scientific name"),
    (111, "Asteroidea", "scientific name"),
    (200, "Placozoa", "scientific name"),
    (201, "Trichoplax", "scientific name"),
    (300, "Streptophyta", "scientific name"),
    (301, "Morus", "scientific name"),
]


def write_taxdump(directory):
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / "nodes.dmp", "w") as handle:
        for taxid, parent, rank in NODES:
            handle.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t-\t|\n")
    with open(directory / "names.dmp", "w") as handle:
        for taxid, name, name_class in NAMES:
            handle.write(f"{taxid}\t|\t{name}\t|\t\t|\t{name_class}\t|\n")
    return directory


class TaxdumpDirMixin:
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.taxdump = str(write_taxdump(Path(cls._tmp.name) / "taxonkit_dbs"))

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()


class CleanSpeciesNameTests(unittest.TestCase):
    """The nominal names carry qualifier noise NCBI will not match."""

    def test_rank_parenthetical_is_stripped(self):
        # Real values from NOVA_260909_LA.
        self.assertEqual(cs.clean_species_name("Asteroidea (Class)"), "Asteroidea")
        self.assertEqual(cs.clean_species_name("Ophiuroidea (Class)"), "Ophiuroidea")

    def test_open_nomenclature_and_backticks_are_stripped(self):
        self.assertEqual(cs.clean_species_name("Astronesthes spp."), "Astronesthes")
        self.assertEqual(
            cs.clean_species_name("Nesogobius sp. `groove cheek`"), "Nesogobius")

    def test_a_clean_binomial_is_left_alone(self):
        self.assertEqual(cs.clean_species_name("Chaceon bicolor"), "Chaceon bicolor")

    def test_empty_input_gives_an_empty_query(self):
        self.assertEqual(cs.clean_species_name(None), "")
        self.assertEqual(cs.clean_species_name("   "), "")


class BuildResolverTests(TaxdumpDirMixin, unittest.TestCase):
    def test_no_taxdump_dir_means_no_resolver(self):
        self.assertIsNone(cs.build_resolver(None))

    def test_a_directory_without_the_dmp_files_means_no_resolver(self):
        # Degrade to species-table-only rather than crashing the run.
        self.assertIsNone(cs.build_resolver("/nonexistent/taxdump"))

    def test_this_pipeline_turns_both_flags_on(self):
        resolver = cs.build_resolver(self.taxdump)
        self.assertTrue(resolver.index_all_ranks)
        self.assertTrue(resolver.class_falls_back_to_phylum)

    def test_an_infraorder_resolves_through_the_built_resolver(self):
        # 'Caridea' is the shape the default five-rank index misses entirely.
        lineage = cs.build_resolver(self.taxdump).lineage_for_name("Caridea")
        self.assertEqual(lineage["class"], "Malacostraca")
        self.assertEqual(lineage["matched_taxid"], 103)


class MatchQueryTests(unittest.TestCase):
    """The lookup reads master_species and matches the cleaned name too."""

    def test_the_maintained_table_is_read_not_species(self):
        self.assertIn("FROM master_species sp", cs.MATCH_QUERY)
        self.assertNotRegex(cs.MATCH_QUERY, r"FROM species\b")

    def test_the_dead_same_genus_fuzzy_tier_is_gone(self):
        # It could never fire: any row it matched, genus_exact had matched first.
        self.assertNotIn("species_fuzzy_same_genus", cs.MATCH_QUERY)

    def test_the_cleaned_name_and_its_genus_are_passed(self):
        params = cs.match_params("Gobiidae spp.")
        self.assertEqual(params, {"raw": "Gobiidae spp.", "clean": "Gobiidae",
                                  "genus": "Gobiidae"})

    def test_a_name_that_cleans_to_nothing_falls_back_to_the_raw_name(self):
        self.assertEqual(cs.match_params("sp.")["clean"], "sp.")

    def test_tiers_are_ordered_species_synonym_then_higher_ranks(self):
        order = [cs.MATCH_QUERY.index(f"'{tier}'") for tier in (
            "species_exact", "synonym_exact", "genus_exact", "family_exact",
            "order_exact", "family_fuzzy", "order_fuzzy")]
        self.assertEqual(order, sorted(order))


class ResolveSpeciesInfoTests(TaxdumpDirMixin, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.resolver = cs.build_resolver(cls.taxdump)

    def resolve(self, db_result, resolver=None):
        with mock.patch.object(cs, "get_species_info", return_value=db_result):
            return cs.resolve_species_info(None, "OG1", resolver)

    def test_an_exact_species_answer_is_left_alone(self):
        result = self.resolve(("Chaceon bicolor", 999, "Curated", "species_exact"),
                              self.resolver)
        self.assertEqual(result,
                         ("Chaceon bicolor", 999, "Curated", "db", "species_exact", ""))

    def test_a_synonym_hit_counts_as_exact(self):
        result = self.resolve(("Old name", 999, "Curated", "synonym_exact"),
                              self.resolver)
        self.assertEqual(result[1:4], (999, "Curated", "db"))

    def test_an_exact_row_without_a_taxid_takes_the_taxdump_one(self):
        # 40% of master_species rows have no ncbi_taxon_id.
        _n, taxon_id, tax_class, source, level, rank = self.resolve(
            ("Chaceon bicolor", None, "Curated", "species_exact"), self.resolver)
        self.assertEqual((taxon_id, tax_class, source), (105, "Curated", "db+taxdump"))
        self.assertEqual((level, rank), ("species_exact", "species"))

    def test_a_congener_row_never_supplies_the_taxid_when_ncbi_has_the_species(self):
        # The bug this replaces: a genus-level hit came back complete, was
        # reported as 'db', and FCS-GX got another species' taxid.
        _n, taxon_id, tax_class, source, level, rank = self.resolve(
            ("Chaceon bicolor", 999, "Curated", "genus_exact"), self.resolver)
        self.assertEqual(taxon_id, 105)
        self.assertEqual(tax_class, "Curated",
                         "the curated class spelling is what BUSCO routing keys on")
        self.assertEqual((source, level, rank), ("db+taxdump", "genus_exact", "species"))

    def test_a_genus_level_hit_takes_the_real_genus_taxid(self):
        _n, taxon_id, tax_class, source, _level, rank = self.resolve(
            ("Chaceon sp 1", 999, "Curated", "genus_exact"), self.resolver)
        self.assertEqual((taxon_id, tax_class, source, rank),
                         (104, "Curated", "db+taxdump", "genus"))

    def test_without_a_taxdump_a_congener_taxid_is_kept_but_flagged(self):
        _n, taxon_id, _class, source, level, rank = self.resolve(
            ("Chaceon bicolor", 999, "Curated", "genus_exact"))
        self.assertEqual((taxon_id, source, level, rank),
                         (999, "db", "genus_exact", ""))

    def test_the_taxdump_beats_a_fuzzy_match(self):
        _n, taxon_id, tax_class, source, level, _rank = self.resolve(
            ("Chaceon", 999, "Wrong", "family_fuzzy"), self.resolver)
        self.assertEqual((taxon_id, tax_class, source), (104, "Malacostraca", "taxdump"))
        self.assertEqual(level, "family_fuzzy")

    def test_a_fuzzy_match_is_the_last_resort(self):
        # 'Chrysogorgidae' for 'Chrysogorgiidae': NCBI has no such name.
        _n, taxon_id, tax_class, source, level, _rank = self.resolve(
            ("Chrysogorgidae", 999, "Anthozoa", "family_fuzzy"), self.resolver)
        self.assertEqual((taxon_id, tax_class, source, level),
                         (999, "Anthozoa", "db", "family_fuzzy"))

    def test_a_missing_row_is_filled_entirely_from_the_taxdump(self):
        nominal, taxon_id, tax_class, source, level, rank = self.resolve(
            ("Chaceon", None, None, None), self.resolver)
        self.assertEqual(nominal, "Chaceon")
        self.assertEqual(taxon_id, 104)
        self.assertEqual(tax_class, "Malacostraca")
        self.assertEqual((source, level, rank), ("taxdump", "", "genus"))

    def test_a_class_name_used_as_the_whole_nominal_name_resolves(self):
        # 'Asteroidea (Class)' is how these arrive from the collection.
        _nominal, taxon_id, tax_class, source, _l, _r = self.resolve(
            ("Asteroidea (Class)", None, None, None), self.resolver)
        self.assertEqual(taxon_id, 111)
        self.assertEqual(tax_class, "Asteroidea")
        self.assertEqual(source, "taxdump")

    def test_a_row_with_no_class_is_topped_up_and_flagged(self):
        _nominal, taxon_id, tax_class, source, _l, _r = self.resolve(
            ("Chaceon", 999, None, "species_exact"), self.resolver)
        self.assertEqual(taxon_id, 999, "an exact row's taxon_id must survive")
        self.assertEqual(tax_class, "Malacostraca")
        self.assertEqual(source, "db+taxdump")

    def test_a_phylum_stands_in_for_a_missing_class(self):
        # BUSCO routes anything unlisted to metazoa anyway.
        _nominal, _taxon_id, tax_class, source, _l, _r = self.resolve(
            ("Trichoplax", None, None, None), self.resolver)
        self.assertEqual(tax_class, "Placozoa")
        self.assertEqual(source, "taxdump")

    def test_a_misspelled_name_stays_unresolved(self):
        # 'Actinaria' for 'Actiniaria' -- no source can rescue a name that does
        # not exist, and guessing one would send FCS-GX the wrong taxon.
        nominal, taxon_id, tax_class, source, _l, _r = self.resolve(
            ("Actinaria", None, None, None), self.resolver)
        self.assertEqual(nominal, "Actinaria")
        self.assertIsNone(taxon_id)
        self.assertIsNone(tax_class)
        self.assertEqual(source, "unresolved")

    def test_a_taxon_outside_metazoa_stays_unresolved(self):
        result = self.resolve(("Morus", None, None, None), self.resolver)
        self.assertEqual(result[3], "unresolved")

    def test_without_a_resolver_an_uncurated_sample_still_fails(self):
        # What a run with no --taxonkit_db_dir must keep getting.
        _nominal, taxon_id, tax_class, source, _l, _r = self.resolve(
            ("Chaceon", None, None, None))
        self.assertIsNone(taxon_id)
        self.assertIsNone(tax_class)
        self.assertEqual(source, "unresolved")

    def test_a_sample_with_no_row_at_all_is_unresolved(self):
        self.assertEqual(self.resolve(None, self.resolver),
                         (None, None, None, "unresolved", "", ""))


class TaxonIdScriptTests(TaxdumpDirMixin, unittest.TestCase):
    """bin/taxon_id.py must answer exactly as the samplesheet does."""

    def run_script(self, db_result):
        sys.path.insert(0, str(BIN))
        try:
            with mock.patch.dict(sys.modules, {"create_samplesheet": cs}):
                spec = importlib.util.spec_from_file_location(
                    "taxon_id_under_test", BIN / "taxon_id.py")
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
        finally:
            sys.path.remove(str(BIN))
        argv = ["taxon_id.py", "/cfg", "OG1", "--taxdump-dir", self.taxdump]
        out = io.StringIO()
        with mock.patch.object(module, "load_db_config", return_value={}), \
             mock.patch.object(cs, "get_species_info", return_value=db_result), \
             mock.patch.object(sys, "argv", argv), \
             contextlib.redirect_stdout(out):
            module.main()
        return out.getvalue().splitlines()

    def test_the_csv_matches_resolve_species_info(self):
        lines = self.run_script(("Chaceon bicolor", 999, "Curated", "genus_exact"))
        self.assertEqual(lines, ["og_id,nominal_species_id,taxon_id,class",
                                 "OG1,Chaceon bicolor,105,Curated"])

    def test_an_unresolved_sample_prints_nothing(self):
        self.assertEqual(self.run_script(("Actinaria", None, None, None)), [])


class ResolutionReportTests(unittest.TestCase):
    def test_every_sample_is_recorded_with_its_source(self):
        # Without this a taxdump-derived class is indistinguishable from a
        # curated one, and there is no way to audit what the fallback did.
        rows = [
            {"sample": "OG1", "nominal_species_id": "Chaceon", "taxon_id": 104,
             "class": "Malacostraca", "source": "taxdump", "match_level": "",
             "taxdump_rank": "genus"},
            {"sample": "OG2", "nominal_species_id": "Actinaria", "taxon_id": "",
             "class": "", "source": "unresolved", "match_level": "",
             "taxdump_rank": ""},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "taxonomy_resolution.tsv"
            cs.write_resolution_report(str(path), rows)
            lines = path.read_text().splitlines()
        self.assertEqual(lines[0].split("\t"), list(cs.RESOLUTION_COLUMNS))
        self.assertEqual(lines[1].split("\t"),
                         ["OG1", "Chaceon", "104", "Malacostraca", "taxdump", "",
                          "genus"])
        self.assertEqual(lines[2].split("\t"),
                         ["OG2", "Actinaria", "", "", "unresolved", "", ""])


class MainWritesTheSheetTests(TaxdumpDirMixin, unittest.TestCase):
    """An unresolved sample must not cost us the whole samplesheet.

    publishDir does not run on a failed task, so exiting non-zero here meant
    `outdir/samplesheet/` stayed empty and there was nothing for an operator to
    inspect or correct. The sheet is now always written, with the bad rows marked
    `unknown`, and the run is stopped by validateTaxonomy() in
    subworkflows/local/prepare_samplesheet instead -- after publishing, before any
    assembly work.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.resolver = cs.build_resolver(cls.taxdump)

    def run_main(self, db_results, tmpdir):
        """Drive main() over a fake set of read pairs, one per db_results key."""
        pairs = [(og, [f"/reads/{og}.ilmn.R1.fq.gz", f"/reads/{og}.ilmn.R2.fq.gz"])
                 for og in db_results]
        argv = ["create_samplesheet.py", "/cfg", "NOVA_260909_LA", "/pooled",
                "--taxdump-dir", self.taxdump,
                "--resolution-report", str(Path(tmpdir) / "taxonomy_resolution.tsv")]
        with mock.patch.object(cs, "discover_pairs_from_dir", return_value=pairs), \
             mock.patch.object(cs, "load_db_config", return_value={}), \
             mock.patch.object(cs, "get_species_info",
                               side_effect=lambda _p, og: db_results[og]), \
             mock.patch.object(sys, "argv", argv):
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                rc = cs.main()
            finally:
                os.chdir(cwd)
        sheet = Path(tmpdir) / "NOVA_260909_LA_samplesheet.csv"
        rows = []
        if sheet.exists():
            with open(sheet, newline="") as handle:
                rows = list(csv.DictReader(handle))
        return rc, rows

    def test_an_unresolved_sample_is_written_as_unknown_and_main_succeeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            rc, rows = self.run_main(
                {"OG1": ("Chaceon", None, None, None),      # rescued by the taxdump
                 "OG2": ("Actinaria", None, None, None)},   # rescued by nothing
                tmp)
        # Exit code is what decides whether publishDir runs at all.
        self.assertIn(rc, (None, 0))
        by_sample = {r["sample"]: r for r in rows}
        self.assertEqual(sorted(by_sample), ["OG1", "OG2"],
                         "the unresolved sample must still get a row")
        self.assertEqual(by_sample["OG2"]["taxon_id"], "unknown")
        self.assertEqual(by_sample["OG2"]["class"], "unknown")
        self.assertEqual(by_sample["OG2"]["nom_species_id"], "Actinaria")

    def test_the_resolved_rows_are_untouched_by_a_bad_neighbour(self):
        with tempfile.TemporaryDirectory() as tmp:
            _rc, rows = self.run_main(
                {"OG1": ("Chaceon", None, None, None),
                 "OG2": ("Actinaria", None, None, None)},
                tmp)
        good = {r["sample"]: r for r in rows}["OG1"]
        self.assertEqual(good["taxon_id"], "104")
        self.assertEqual(good["class"], "Malacostraca")

    def test_the_placeholder_matches_what_the_input_schema_allows(self):
        # schema_input.json permits taxon_id '^(unknown|None)$', so these rows
        # parse and reach validateTaxonomy rather than failing validation first.
        schema = json.load(open(ROOT / "assets" / "schema_input.json"))
        pattern = [b["pattern"] for b in
                   schema["items"]["properties"]["taxon_id"]["anyOf"] if "pattern" in b]
        self.assertTrue(re.match(pattern[0], cs.UNKNOWN))

    def test_no_usable_read_pairs_still_fails(self):
        # A different failure from unresolved taxonomy, and there is no sheet to
        # write, so this one must keep its non-zero exit.
        with tempfile.TemporaryDirectory() as tmp:
            rc, rows = self.run_main({}, tmp)
        self.assertEqual(rc, 1)
        self.assertEqual(rows, [])


class BuscoVertebrateClassTests(unittest.TestCase):
    """Every vertebrate class master_species writes must reach the vertebrata lineage.

    The class comes from master_species whenever it matches, and it spells sharks
    'Elasmobranchii', not NCBI's 'Chondrichthyes'. While buscoDbForClass only knew
    the NCBI name, every shark was scored against metazoa.
    """

    def test_master_species_spellings_are_routed_to_vertebrata(self):
        source = (ROOT / "subworkflows" / "local" /
                  "utils_nfcore_oceangenomes_draftgenomes_pipeline" / "main.nf").read_text()
        block = re.search(r"BUSCO_VERT_CLASSES = \[(.*?)\]", source, re.S).group(1)
        routed = set(re.findall(r"'([^']+)'", block))
        for name in ("Elasmobranchii", "Holocephali", "Petromyzonti", "Coelacanthi",
                     "Crocodylia", "Myxini", "Mammalia", "Aves", "Reptilia"):
            self.assertIn(name, routed)


class SchemaMetaAnnotationTests(unittest.TestCase):
    """The samplesheet is read by column name, not by position.

    prepare_samplesheet used to rebuild its meta map from sample_record[4], [5]
    and [6]; reordering a column in schema_input.json would have silently
    reassigned every taxonomy field.
    """

    def setUp(self):
        self.schema = json.load(open(ROOT / "assets" / "schema_input.json"))["items"]

    def test_every_meta_column_is_annotated(self):
        props = self.schema["properties"]
        expected = {"sample": "id", "run": "run", "date": "date", "prefix": "prefix",
                    "nom_species_id": "nom_species_id", "taxon_id": "taxon_id",
                    "class": "class"}
        for column, key in expected.items():
            self.assertEqual(props[column].get("meta"), [key],
                             f"{column} must be annotated as meta {key}")

    def test_the_read_columns_are_not_meta(self):
        # They are the tuple payload; folding them into meta would change its shape.
        for column in ("fastq_1", "fastq_2"):
            self.assertIsNone(self.schema["properties"][column].get("meta"))

    def test_every_column_stays_required(self):
        # This is what makes the meta annotation safe. nf-schema substitutes an
        # EMPTY LIST for an absent optional column, which changes the shape of
        # every meta map and silently defeats join(..., by: 0) on resume. With no
        # optional columns that placeholder can never appear.
        self.assertEqual(sorted(self.schema["required"]),
                         sorted(self.schema["properties"]))


if __name__ == "__main__":
    unittest.main()
