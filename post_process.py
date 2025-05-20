#!/usr/bin/env python3
# /// script
# requires-python = ">=3.7"
# dependencies = [
#   "pandas",
#   "tqdm",
# ]
# ///

# SPECTRONAUT POST PROCESSING
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


# version tracking
__version = "1.1.5"
__date = "2025-03-26"

# PARAMETERS

# these parameters should be set accordingly
CROSSLINKER = "PhoX" # name of the crosslinker
CROSSLINKER_MASS = 209.97181 # delta mass of the crosslinker
SPECTRONAUT_DELIM = "," # delimiter in Spectronaut output file, e.g. "," for comma delimited files, "\t" for tab delimited files
SPECTRONAUT_MATCH_TOLERANCE = 0.05 # match tolerance in Da
SPECTRONAUT_FRAGMENT_MZ_COLUMN_NAME = "F.CalibratedMz" # which F Mz to use for matching
SPECTRONAUT_CSCORE_COLUMN_NAME = "EG.Cscore" # which Cscore to use for re-soring

# import packages
import argparse
import pandas as pd
from tqdm import tqdm

from typing import Any
from typing import Dict
from typing import List

##################### POST PROCESS #####################

def get_mz_key(mz: float) -> float:
    #return str(int(mz * 1000))
    return mz

def get_fragment_key(mz: float) -> str:
    return f"{round(mz, 4):.4f}"

def get_key_spec_lib(row: pd.Series) -> str:
    # ModifiedPeptide
    # DAKQRIVDK_NGVKM[Oxidation]C[Carbamidomethyl]PR
    # PrecursorCharge
    # 4
    # linkId
    # A0A4V0IIP8_Q21298-126_45
    return f"{str(row['ModifiedPeptide']).strip()}.{int(row['PrecursorCharge'])}:{str(row['linkId']).strip()}"

def get_key_spectronaut(row: pd.Series) -> str:
    # EG.PrecursorId
    # AAHHADGLAKGLHETC[Carbamidomethyl]K_M[Oxidation]FIPKSHTK.5
    # PG.ProteinNames
    # P10771_P10771-113_46
    # ---> if PG.ProteinNames
    # P10771_P10771
    # use FG.Comment
    # 113_46
    if "-" in str(row["PG.ProteinNames"]):
        return f"{str(row['EG.PrecursorId']).strip()}:{str(row['PG.ProteinNames']).strip()}"
    return f"{str(row['EG.PrecursorId']).strip()}:{str(row['PG.ProteinNames']).strip()}-{str(row['FG.Comment']).strip()}"

def read_spectral_library(filename: str) -> Dict[str, Dict[str, Any]]:

    spec_lib = pd.read_csv(filename, low_memory = False)

    index = dict()
    for i, row in tqdm(spec_lib.iterrows(), total = spec_lib.shape[0], desc = "Creating index..."):
        key = get_key_spec_lib(row)
        if key in index:
            index[key]["rows"].append(row)

            if int(row["FragmentPepId"]) == 0:
                index[key]["total_ions_a"] += 1
            else:
                index[key]["total_ions_b"] += 1

            ion_mz = get_mz_key(float(row["FragmentMz"]))
            if ion_mz in index[key]["ions"]:
                index[key]["ions"][ion_mz].append(row)
            else:
                index[key]["ions"][ion_mz] = [row]
        else:
            index[key] = {"rows": [row],
                          "total_ions_a": 1 if int(row["FragmentPepId"]) == 0 else 0,
                          "total_ions_b": 1 if int(row["FragmentPepId"]) == 1 else 0,
                          "ions": {get_mz_key(float(row["FragmentMz"])): [row]}}

    return index

def generate_fragment_index(spectronaut: pd.DataFrame, index: dict) -> Dict[str, Dict[str, Any]]:
    # dict keeping track of annotated fragments for every unique crosslink id
    fragment_annotation = dict()
    for i, row in tqdm(spectronaut.iterrows(), total = spectronaut.shape[0], desc = "Generating fragment ion index..."):
        # unique crosslink id
        key = get_key_spectronaut(row)
        # current fragment ion from spectronaut row
        ion = float(row[SPECTRONAUT_FRAGMENT_MZ_COLUMN_NAME])
        # all spectral library ions for the crosslink id
        # this is a dict that matches ion mass -> list(fragments_with_that_mass[pd.Series])
        ions = index[key]["ions"]
        # for every ion mass of the crosslink id in the spec lib, do:
        for current_ion_mz in ions.keys():
            fragment_key_part = get_fragment_key(current_ion_mz)
            # check if spectronaut fragment is within tolerance bounds of spec lib ion
            if round(current_ion_mz, 4) > round(ion - SPECTRONAUT_MATCH_TOLERANCE, 4) and round(current_ion_mz, 4) < round(ion + SPECTRONAUT_MATCH_TOLERANCE, 4):
                # if fragment annotation struct does not have entry for crosslink id, create one that is empty
                if key not in fragment_annotation:
                    fragment_annotation[key] = {"matched_number_ions_a": 0,
                                                "matched_number_ions_b": 0,
                                                "fragments": []}
                # for all ions that match mass, do:
                for current_ion in ions[current_ion_mz]:
                    # if ion of peptide a
                    if current_ion["FragmentPepId"] == 0:
                        # generate a unique fragment key
                        fragment_key = fragment_key_part + "_0"
                        # check if a fragment of this type has not already been annotated, we don't want to annotate the same ion twice
                        if fragment_key not in fragment_annotation[key]["fragments"]:
                            fragment_annotation[key]["matched_number_ions_a"] += 1
                            fragment_annotation[key]["fragments"].append(fragment_key)
                    else:
                        # same for peptide b
                        fragment_key = fragment_key_part + "_1"
                        if fragment_key not in fragment_annotation[key]["fragments"]:
                            fragment_annotation[key]["matched_number_ions_b"] += 1
                            fragment_annotation[key]["fragments"].append(fragment_key)
            # we do not break here, because we want to check the rest of the spec lib ions in case the spectronaut ion matches
            # for both peptide a and peptide b

    return fragment_annotation

def annotate_spectronaut_result(filename: str) -> pd.DataFrame:

    spectronaut = pd.read_csv(filename, sep = SPECTRONAUT_DELIM, low_memory = False)
    filename_spec_lib = str(spectronaut["EG.Library"].at[0])
    index = read_spectral_library(filename_spec_lib)

    ## available columns in spec lib
    # COLUMN                                            EXAMPLE
    # linkId                                            P18948_P18948-516_509
    # ProteinID                                         P18948_P18948
    # Organism                                          Caenorhabditis elegans
    # StrippedPeptide                                   DAKQRIVDKNGVKMCPR
    # FragmentGroupId                                   DAKQRIVDK_NGVKMCPR-3_4:5
    # PrecursorCharge                                   5
    # PrecursorMz                                       429.823708
    # ModifiedPeptide                                   DAKQRIVDK_NGVKM[Oxidation]C[Carbamidomethyl]PR
    # IsotopeLabel                                      0
    # File                                              20240131_E0_NEO6_Mueller_MS_TechHub_IMP_THIDVJSSG004_Celegans_Nulcei_S3_DSG_SEC9_allng_FAIMS_001.raw
    # scanID                                            15497
    # run                                               20240131_E0_NEO6_THIDVJSSG004_Celegans_Nulcei_S1toS3_DSG_SEC_entrapment
    # searchID                                          MS Annika
    # crosslinkedResidues                               516_509
    # LabeledSequence                                   DAKQRIVDK_NGVKM[Oxidation]C[Carbamidomethyl]PR
    # iRT                                               -17.5246511627907
    # RT                                                32.1844
    # CCS                                               0
    # IonMobility                                       -50
    # FragmentCharge                                    1
    # FragmentType                                      b
    # FragmentNumber                                    2
    # FragmentPepId                                     0
    # FragmentMz                                        187.0715942
    # RelativeIntensity                                 0.115658057537211
    # FragmentLossType                                  None
    # CLContainingFragment                              False
    # LossyFragment                                     False
    # isDecoy                                           False
    # DecoyType                                         TT
    # PeptidePositions                                  503_497

    # for rescoring
    # a: matched number ions a
    # b: total number ions a
    # c: matched number ions b
    # d: total number ions b
    # e: relative match score a (a / b)
    # f: relative match score b (c / d)
    # g: partial c score a (c score * (a / a + c))
    # h: partial c score b (c score * (c / a + c))
    # i: composite relative match score min(e, f)
    # j: composite partial c score min(g, h)
    fragment_annotation = generate_fragment_index(spectronaut, index)

    # a
    def annotate_MatchedIonsA(row: pd.Series, fragment_annotation: dict) -> int:
        key = get_key_spectronaut(row)
        return fragment_annotation[key]["matched_number_ions_a"]

    tqdm.pandas(desc = "Annotating number of matched ions A...")
    spectronaut["PP.MatchedIonsA"] = spectronaut.progress_apply(lambda row: annotate_MatchedIonsA(row, fragment_annotation), axis = 1)

    # b
    def annotate_TotalIonsA(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return index[key]["total_ions_a"]

    tqdm.pandas(desc = "Annotating number of total ions A...")
    spectronaut["PP.TotalIonsA"] = spectronaut.progress_apply(lambda row: annotate_TotalIonsA(row, index), axis = 1)

    # c
    def annotate_MatchedIonsB(row: pd.Series, fragment_annotation: dict) -> int:
        key = get_key_spectronaut(row)
        return fragment_annotation[key]["matched_number_ions_b"]

    tqdm.pandas(desc = "Annotating number of matched ions B...")
    spectronaut["PP.MatchedIonsB"] = spectronaut.progress_apply(lambda row: annotate_MatchedIonsB(row, fragment_annotation), axis = 1)

    # d
    def annotate_TotalIonsB(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return index[key]["total_ions_b"]

    tqdm.pandas(desc = "Annotating number of total ions B...")
    spectronaut["PP.TotalIonsB"] = spectronaut.progress_apply(lambda row: annotate_TotalIonsB(row, index), axis = 1)

    # e
    tqdm.pandas(desc = "Annotating relative match score A...")
    spectronaut["PP.RelativeMatchScoreA"] = spectronaut.progress_apply(lambda row: row["PP.MatchedIonsA"] / row["PP.TotalIonsA"], axis = 1)

    # f
    tqdm.pandas(desc = "Annotating relative match score B...")
    spectronaut["PP.RelativeMatchScoreB"] = spectronaut.progress_apply(lambda row: row["PP.MatchedIonsB"] / row["PP.TotalIonsB"], axis = 1)

    # g
    def annotate_PartialCscoreA(row: pd.Series) -> float:
        cscore = row[SPECTRONAUT_CSCORE_COLUMN_NAME]
        partial = row["PP.MatchedIonsA"] / (row["PP.MatchedIonsA"] + row["PP.MatchedIonsB"])
        return cscore * partial

    tqdm.pandas(desc = "Annotating partial Cscore A...")
    spectronaut["PP.PartialCscoreA"] = spectronaut.progress_apply(lambda row: annotate_PartialCscoreA(row), axis = 1)

    # h
    def annotate_PartialCscoreB(row: pd.Series) -> float:
        cscore = row[SPECTRONAUT_CSCORE_COLUMN_NAME]
        partial = row["PP.MatchedIonsB"] / (row["PP.MatchedIonsA"] + row["PP.MatchedIonsB"])
        return cscore * partial

    tqdm.pandas(desc = "Annotating partial Cscore B...")
    spectronaut["PP.PartialCscoreB"] = spectronaut.progress_apply(lambda row: annotate_PartialCscoreB(row), axis = 1)

    # i
    tqdm.pandas(desc = "Annotating composite relative match score...")
    spectronaut["PP.CompositeRelativeMatchScore"] = spectronaut.progress_apply(lambda row: min(row["PP.RelativeMatchScoreA"], row["PP.RelativeMatchScoreB"]), axis = 1)

    # j
    tqdm.pandas(desc = "Annotating composite partial Cscore...")
    spectronaut["PP.CompositePartialCscore"] = spectronaut.progress_apply(lambda row: min(row["PP.PartialCscoreA"], row["PP.PartialCscoreB"]), axis = 1)

    # annotation of other properties
    def annotate_DecoyType(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key]["rows"][0]["DecoyType"]).strip()

    tqdm.pandas(desc = "Annotating decoy type...")
    spectronaut["PP.DecoyType"] = spectronaut.progress_apply(lambda row: annotate_DecoyType(row, index), axis = 1)

    def annotate_ProteinA(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key]["rows"][0]["ProteinID"].split("_")[0]).strip()

    tqdm.pandas(desc = "Annotating protein A...")
    spectronaut["PP.ProteinA"] = spectronaut.progress_apply(lambda row: annotate_ProteinA(row, index), axis = 1)

    def annotate_ProteinB(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key]["rows"][0]["ProteinID"].split("_")[1]).strip()

    tqdm.pandas(desc = "Annotating protein B...")
    spectronaut["PP.ProteinB"] = spectronaut.progress_apply(lambda row: annotate_ProteinB(row, index), axis = 1)

    def annotate_CrosslinkPositionProteinA(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return int(index[key]["rows"][0]["linkId"].split("-")[1].split("_")[0])

    tqdm.pandas(desc = "Annotating crosslink position in protein A...")
    spectronaut["PP.CrosslinkPositionProteinA"] = spectronaut.progress_apply(lambda row: annotate_CrosslinkPositionProteinA(row, index), axis = 1)

    def annotate_CrosslinkPositionProteinB(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return int(index[key]["rows"][0]["linkId"].split("-")[1].split("_")[1])

    tqdm.pandas(desc = "Annotating crosslink position in protein B...")
    spectronaut["PP.CrosslinkPositionProteinB"] = spectronaut.progress_apply(lambda row: annotate_CrosslinkPositionProteinB(row, index), axis = 1)

    def annotate_PeptideA(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key]["rows"][0]["FragmentGroupId"].split("-")[0].split("_")[0]).strip()

    tqdm.pandas(desc = "Annotating peptide sequence A...")
    spectronaut["PP.PeptideA"] = spectronaut.progress_apply(lambda row: annotate_PeptideA(row, index), axis = 1)

    def annotate_PeptideB(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key]["rows"][0]["FragmentGroupId"].split("-")[0].split("_")[1]).strip()

    tqdm.pandas(desc = "Annotating peptide sequence B...")
    spectronaut["PP.PeptideB"] = spectronaut.progress_apply(lambda row: annotate_PeptideB(row, index), axis = 1)

    def annotate_CrosslinkPositionPeptideA(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return int(index[key]["rows"][0]["FragmentGroupId"].split("-")[1].split(":")[0].split("_")[0])

    tqdm.pandas(desc = "Annotating crosslink position in peptide A...")
    spectronaut["PP.CrosslinkPositionPeptideA"] = spectronaut.progress_apply(lambda row: annotate_CrosslinkPositionPeptideA(row, index), axis = 1)

    def annotate_CrosslinkPositionPeptideB(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return int(index[key]["rows"][0]["FragmentGroupId"].split("-")[1].split(":")[0].split("_")[1])

    tqdm.pandas(desc = "Annotating crosslink position in peptide B...")
    spectronaut["PP.CrosslinkPositionPeptideB"] = spectronaut.progress_apply(lambda row: annotate_CrosslinkPositionPeptideB(row, index), axis = 1)

    def annotate_PeptidoformA(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key]["rows"][0]["ModifiedPeptide"].split("_")[0]).strip()

    tqdm.pandas(desc = "Annotating peptidoform sequence A...")
    spectronaut["PP.PeptidoformA"] = spectronaut.progress_apply(lambda row: annotate_PeptidoformA(row, index), axis = 1)

    def annotate_PeptidoformB(row: pd.Series, index: dict) -> str:
        key = get_key_spectronaut(row)
        return str(index[key]["rows"][0]["ModifiedPeptide"].split("_")[1]).strip()

    tqdm.pandas(desc = "Annotating peptidoform sequence B...")
    spectronaut["PP.PeptidoformB"] = spectronaut.progress_apply(lambda row: annotate_PeptidoformB(row, index), axis = 1)

    def annotate_PeptidePositionA(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return int(index[key]["rows"][0]["PeptidePositions"].split("_")[0])

    tqdm.pandas(desc = "Annotating peptide A position in protein A...")
    spectronaut["PP.PeptidePositionProteinA"] = spectronaut.progress_apply(lambda row: annotate_PeptidePositionA(row, index), axis = 1)

    def annotate_PeptidePositionB(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return int(index[key]["rows"][0]["PeptidePositions"].split("_")[1])

    tqdm.pandas(desc = "Annotating peptide B position in protein B...")
    spectronaut["PP.PeptidePositionProteinB"] = spectronaut.progress_apply(lambda row: annotate_PeptidePositionB(row, index), axis = 1)

    def annotate_DecoyA(row: pd.Series, index: dict) -> bool:
        key = get_key_spectronaut(row)
        d = str(index[key]["rows"][0]["DecoyType"]).strip()[0]
        if d not in ["T", "D"]:
            raise ValueError(f"DecoyType not recognized: {d}")
        return d == "D"

    tqdm.pandas(desc = "Annotating decoy type for peptide A...")
    spectronaut["PP.IsDecoyA"] = spectronaut.progress_apply(lambda row: annotate_DecoyA(row, index), axis = 1)

    def annotate_DecoyB(row: pd.Series, index: dict) -> bool:
        key = get_key_spectronaut(row)
        d = str(index[key]["rows"][0]["DecoyType"]).strip()[1]
        if d not in ["T", "D"]:
            raise ValueError(f"DecoyType not recognized: {d}")
        return d == "D"

    tqdm.pandas(desc = "Annotating decoy type for peptide B...")
    spectronaut["PP.IsDecoyB"] = spectronaut.progress_apply(lambda row: annotate_DecoyB(row, index), axis = 1)

    def annotate_SourceScanID(row: pd.Series, index: dict) -> int:
        key = get_key_spectronaut(row)
        return int(index[key]["rows"][0]["scanID"])

    tqdm.pandas(desc = "Annotating spectral library source scan ID...")
    spectronaut["PP.SourceScanID"] = spectronaut.progress_apply(lambda row: annotate_SourceScanID(row, index), axis = 1)

    spectronaut["PP.PseudoScanNumber"] = pd.Series(range(spectronaut.shape[0]))
    spectronaut["PP.Crosslinker"] = pd.Series([CROSSLINKER for i in range(spectronaut.shape[0])])
    spectronaut["PP.CrosslinkerMass"] = pd.Series([CROSSLINKER_MASS for i in range(spectronaut.shape[0])])

    return spectronaut

def group_by_residue_pair(data: pd.DataFrame) -> pd.DataFrame:
    rows = list()
    seen_residue_pairs = set()

    for i, row in tqdm(data.iterrows(), total = data.shape[0], desc = "Grouping results by residue pairs..."):
        key = get_key_spectronaut(row)
        if key not in seen_residue_pairs:
            seen_residue_pairs.add(key)
            rows.append(row)

    return pd.concat(rows, axis = 1).T

def export_to_xiFDR(data: pd.DataFrame) -> pd.DataFrame:
    # col matching
    # run <- R.Condition, R.FileName
    # scan <- PP.PseudoScanNumber
    # peptide1 <- PP.PeptideA
    # peptide2 <- PP.PeptideB
    # peptide link 1 <- PP.CrosslinkPositionPeptideA
    # peptide link 2 <- PP.CrosslinkPositionPeptideB
    # is decoy 1 <- PP.IsDecoyA (-> "true", "false")
    # is decoy 2 <- PP.IsDecoyB (-> "true", "false")
    # precursor charge <- FG.Charge
    # accession1 <- PP.ProteinA
    # accession2 <- PP.ProteinB
    # peptide position 1 <- PP.PeptidePositionProteinA
    # peptide position 2 <- PP.PeptidePositionProteinB
    # score <- PP.CompositeRelativeMatchScore || PP.CompositePartialCscore

    data["run"] = data.apply(lambda row: f"{row['R.Condition'].strip()}:{row['R.FileName'].strip()}", axis = 1)
    data.rename(columns = {"PP.PseudoScanNumber": "scan"}, inplace = True)
    data.rename(columns = {"PP.PeptideA": "peptide1"}, inplace = True)
    data.rename(columns = {"PP.PeptideB": "peptide2"}, inplace = True)
    data.rename(columns = {"PP.CrosslinkPositionPeptideA": "peptide link 1"}, inplace = True)
    data.rename(columns = {"PP.CrosslinkPositionPeptideB": "peptide link 2"}, inplace = True)
    data["is decoy 1"] = data.apply(lambda row: "true" if bool(row["PP.IsDecoyA"]) else "false", axis = 1)
    data["is decoy 2"] = data.apply(lambda row: "true" if bool(row["PP.IsDecoyB"]) else "false", axis = 1)
    data["precursor charge"] = data.apply(lambda row: int(row["FG.Charge"]), axis = 1)
    data.rename(columns = {"PP.ProteinA": "accession1"}, inplace = True)
    data.rename(columns = {"PP.ProteinB": "accession2"}, inplace = True)
    data.rename(columns = {"PP.PeptidePositionProteinA": "peptide position 1"}, inplace = True)
    data.rename(columns = {"PP.PeptidePositionProteinB": "peptide position 2"}, inplace = True)
    data.rename(columns = {"PP.Crosslinker": "Crosslinker"}, inplace = True)
    data.rename(columns = {"PP.CrosslinkerMass": "CrosslinkerMass"}, inplace = True)

    cols_to_keep = ["run", "scan", "peptide1", "peptide2", "peptide link 1", "peptide link 2",
                    "is decoy 1", "is decoy 2", "precursor charge", "accession1", "accession2",
                    "peptide position 1", "peptide position 2", "PP.CompositeRelativeMatchScore", "PP.CompositePartialCscore",
                    "Crosslinker", "CrosslinkerMass",
                    SPECTRONAUT_CSCORE_COLUMN_NAME]
    return data[cols_to_keep]

def main(argv = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(metavar = "f",
                        dest = "file",
                        help = "Name/Path of the Spectronaut result file to process.",
                        type = str,
                        nargs = 1)
    parser.add_argument("-y", "--non-interactive",
                        dest = "non_interactive",
                        action = "store_true",
                        default = False,
                        help = "Skip information check.")
    parser.add_argument("--version",
                        action = "version",
                        version = __version)
    args = parser.parse_args(argv)

    filename = args.file[0]

    if args.non_interactive:
        r = annotate_spectronaut_result(filename)
    else:
        y = input("Some Spectronaut specific parameters have to be set in the python script. " +
                  "Please refer to the documentation and confirm that you set the parameters accordingly [y/n]:")
        if y.lower().strip() not in ["y", "yes"]:
            raise RuntimeError("Post processing was terminated because of missing confirmation!")
        r = annotate_spectronaut_result(filename)

    output_1 = filename + "_annotated.csv"
    output_2 = output_1 + "_grouped_by_residue_pair.csv"
    output_3 = output_2 + "_xiFDR.csv"

    print("Writing annotated Spectronaut result to file...")
    r.to_csv(output_1, index = False)
    print(f"Finished writing {output_1}.")
    print("Writing grouped Spectronaut result to file...")
    g = group_by_residue_pair(r)
    g.to_csv(output_2, index = False)
    print(f"Finished writing {output_2}.")
    print("Writing grouped Spectronaut result in xiFDR format to file...")
    x = export_to_xiFDR(g)
    x.to_csv(output_3, index = False)
    print(f"Finished writing {output_3}.")
    print("Finished post processing!")

    return 0

if __name__ == "__main__":

    print(main())
