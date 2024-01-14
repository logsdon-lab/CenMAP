# This script is used after HMMER performed alpha satellite HOR analysis
# Convert HMMER table output to BED format
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
#
# Usage: awk -v chr_name=chrN -v coords=[ali|env] -v th=0.8 -v skipsf=1
# -f hmmertblout2bed.awk input-hmmer-tbl.out > output.bed

BEGIN {
  OFS = "\t";
  IGNORECASE = 1;
  if (coords == "env") coordsEnv = 1;
  if (+th < 0) th = 0;
  if (+th > 1) th = 1;

  cnA["M1"]                    = "255,255,0";
  cnA["Qa"]                    = "165,39,0";
  cnA["Pa"]                    = "75,59,162";
  cnA["Ta"]                    = "75,59,162";
  cnA["Ea"]                    = "153,153,255";
  cnA["Fa"]                    = "0,128,128";
  cnA["Ia"]                    = "153,51,102";
  cnA["Aa"]                    = "172,172,172";
  cnA["Ja"]                    = "225,126,231";
  cnA["Ba"]                    = "255,153,0";
  cnA["Ca"]                    = "224,0,64";
  cnA["Oa"]                    = "255,229,153";
  cnA["Na"]                    = "32,160,64";
  cnA["Ka"]                    = "0,255,0";
  cnA["Ha"]                    = "127,96,0";
  cnA["Ga"]                    = "255,255,0";
  cnA["R1"]                    = "0,96,192";
  cnA["R2"]                    = "102,153,255";
  cnA["D1"]                    = "153,0,255";
  cnA["D2"]                    = "210,110,250";
  cnA["D3"]                    = "255,187,255";
  cnA["D4"]                    = "205,150,205";
  cnA["D5"]                    = "181,102,195";
  cnA["D6"]                    = "159,84,183";
  cnA["D7"]                    = "176,63,176";
  cnA["D8"]                    = "141,18,141";
  cnA["D9"]                    = "102,30,102";
  cnA["FD"]                    = "139,102,139";
  cnA["W"]                     = "0,255,255";
  cnA["J1"]                    = "234,153,153";
  cnA["J2"]                    = "255,204,204";
  cnA["J3"]                    = "255,20,147";
  cnA["J4"]                    = "205,16,118";
  cnA["J5"]                    = "205,96,144";
  cnA["J6"]                    = "205,140,149";
  cnA["La"]                    = "128,0,128";
  cnA["S4/6C13/14/21/22H2\\."] = "255,10,10";
  cnA["S4/6C13/14/21H1\\."]    = "255,10,10";
  cnA["S4/6C13/14/21/22H8\\."] = "99,171,63";
  cnA["S4C13/14/21/22H4\\."]   = "94,255,148";
  cnA["S4C13/14/21/22H5\\."]   = "7,171,171";
  cnA["S4C22H3\\."]            = "7,171,171";
  cnA["S4C13/14/21/22H9\\."]   = "10,92,255";
  cnA["S4C15H2\\."]            = "120,120,171";
  cnA["S4C15H2-B\\."]          = "230,179,255";
  cnA["S4C15H3\\."]            = "230,179,255";
  cnA["S4C15H3-A\\."]          = "230,179,255";
  cnA["S4C15H3-B\\."]          = "171,120,137";
  cnA["S4C15H3-AB\\."]         = "180,180,180";
  cnA["S4C20H4\\."]            = "171,63,135";
#  cnA["S4C20H5-A\\."]          = "255,173,10";
  cnA["S4C20H7\\."]            = "255,173,10";
#  cnA["S4C20H5-B\\."]          = "94,255,255";
  cnA["S4C20H8\\."]            = "94,255,255";
  cnA["S4C20H5-C\\."]          = "7,62,171";
  cnA["S4C20H5-D\\."]          = "92,10,255";
  cnA["S4CYH1L\\."]            = "171,7,7";
  cnA["S6C13/14/21/22H3\\."]   = "255,179,230";
  cnA["S6C13/14/21/22H3-A\\."] = "255,179,230";
  cnA["S6C22H2-A\\."]          = "255,179,230";
  cnA["S6C13/14/21/22H3-B\\."] = "255,94,94";
  cnA["S6C22H2-B\\."]          = "255,94,94";
  cnA["S5C11H3\\."]            = "148,94,255";
  cnA["S5C11H4\\."]            = "148,94,255";
  cnA["S5C5pH5\\."]            = "171,116,7";
  cnA["S5C5H5\\."]             = "171,116,7";
  cnA["S5C5/19pH5\\."]         = "171,116,7";
  cnA["S5C5/19H5\\."]          = "171,116,7";
  cnA["S5C5pH6\\."]            = "173,255,10";
  cnA["S5C5H6\\."]             = "173,255,10";
  cnA["S5C5pH7-A\\."]          = "137,171,120";
  cnA["S5C5H7-A\\."]           = "137,171,120";
  cnA["S5C5/19pH7-A\\."]       = "137,171,120";
  cnA["S5C5/19H7-A\\."]        = "137,171,120";
  cnA["S5C5pH7-B\\."]          = "255,179,204";
  cnA["S5C5H7-B\\."]           = "255,179,204";
  cnA["S5C5/19pH7-B\\."]       = "255,179,204";
  cnA["S5C5/19H7-B\\."]        = "255,179,204";
  cnA["S5C7H2\\."]             = "178,255,204";
  cnA["S5C5/19qH4-A\\."]       = "63,171,171";
  cnA["S5C5/19H4-A\\."]        = "63,171,171";
  cnA["S5C5/19qH4-B\\."]       = "255,10,255";
  cnA["S5C5/19H4-B\\."]        = "255,10,255";
  cnA["S5C13/14/21/22H6\\."]   = "171,63,63";
  cnA["S5C13/14/21H2\\."]      = "171,63,63";
  cnA["S5C20H6\\."]            = "255,201,94";
  cnA["S5C4H2\\."]             = "135,63,171"; #"116,171,7";
  cnA["S5C5/19H9d\\."]         = "255,92,10";
  cnA["S5C1qH4d\\."]           = "255,10,92";
  cnA["S5C1H4d\\."]            = "255,10,92";
  cnA["S5C1qH5d\\."]           = "10,255,10";
  cnA["S5C1H5d\\."]            = "10,255,10";
  cnA["S5C1qH6d\\."]           = "63,99,171";
  cnA["S5C1H6d\\."]            = "63,99,171";
  cnA["S3C1pH2-B\\."]          = "171,7,171";
  cnA["S3C1H2-B\\."]           = "171,7,171";
  cnA["S3C1pH2-A\\."]          = "62,7,171"; # "120,171,137";
  cnA["S3C1H2-A\\."]           = "62,7,171"; # "120,171,137";
  cnA["S3C1qH2-D\\."]          = "255,179,179";
  cnA["S3C1H2-D\\."]           = "255,179,179";
  cnA["S3C1qH2-C\\."]          = "171,135,63";
  cnA["S3C1H2-C\\."]           = "171,135,63";
  cnA["S3C11H1L\\."]           = "201,255,94";
  # cnA["S3C11H2\\."]            = "7,171,7";      # FEDOR: renamed to S3C11H2-A
  cnA["S3CXH1L\\."]            = "10,255,173";
  cnA["S3C17H1-B\\."]          = "120,171,171";
  cnA["S3C17H1L\\."]           = "179,204,255";
  cnA["S3C17H1-C\\."]          = "99,63,171";
  cnA["S3C17H2\\."]            = "255,94,255";
  cnA["S3C17H2d\\."]           = "255,94,255";
  cnA["S3C11H3d\\."]           = "171,7,62";
  cnA["S2C2H1L\\."]            = "171,120,120";
  cna["S2C2qH2-A\\."]          = "255,148,94";
  cna["S2C2H2-A\\."]           = "255,148,94";
  cnA["S2C2pH2-B\\."]          = "255,230,179";
  cnA["S2C2H2-B\\."]           = "255,230,179";
  cnA["S2C4H1L\\."]            = "135,171,63";
  cnA["S2C8H1L\\."]            = "94,255,94";
  cnA["S2C9H1L\\."]            = "7,171,116";
  cnA["S2C15H1L\\."]           = "10,173,255";
  cnA["S2C13/21H1L\\."]        = "120,137,171";
  cnA["S2C13/21H1-B\\."]       = "204,179,255";
  cnA["S2C14/22H1L\\."]        = "171,63,171";
  cnA["S2C16pH2-B\\."]         = "63,171,63";
  cnA["S2C16H2-B\\."]          = "63,171,63";
  cnA["S2C16pH2-A\\."]         = "255,94,148";
  cnA["S2C16H2-A\\."]          = "255,94,148";
  cnA["S2C16pH3d\\."]          = "171,154,120";
  cnA["S2C16H3d\\."]           = "171,154,120";
  cnA["S2C16pH4d\\."]          = "94,255,201";
  cnA["S2C16H4d\\."]           = "94,255,201";
  cnA["S2CMH4d\\."]            = "94,255,201";
  cnA["S2C18pH2-A\\."]         = "7,116,171";
  cnA["S2C18H2-A\\."]          = "7,116,171";
  cnA["S2C18H1L\\."]           = "10,10,255";
  cnA["S2C18qH2-B\\."]         = "255,179,255";
  cnA["S2C18H2-B\\."]          = "255,179,255";
  cnA["S2C18qH2-D\\."]         = "171,63,99";
  cnA["S2C18H2-D\\."]          = "171,63,99";
  cnA["S2C18qH2-C\\."]         = "94,201,255";
  cnA["S2C18H2-C\\."]          = "94,201,255";
  cnA["S2C18qH2-E\\."]         = "173,10,255";
  cnA["S2C18H2-E\\."]          = "173,10,255";
  cnA["S2C20H2\\."]            = "154,171,120";
  cnA["S2C20H1L\\."]           = "7,7,171";
  cnA["S2C20H3\\."]            = "178,255,178";
  cnA["S02C20H3\\."]           = "178,255,178";
  cnA["S2CMH2d\\."]            = "63,171,135";
  cnA["S02CMH2d\\."]           = "63,171,135";
  cnA["S2C20H5d\\."]           = "171,120,171";
  cnA["S1C3H1L\\."]            = "171,171,7";
  cnA["S01/1C3H1L\\."]         = "171,171,7";
  cnA["S1C3H2\\."]             = "92,255,10";
  cnA["S01C3H2\\."]            = "92,255,10";
  cnA["S1C3H3d\\."]            = "63,63,171";
  cnA["S01C3/6H1d\\."]         = "63,63,171";
  cnA["S1CMH1d\\."]            = "201,94,255";
  cnA["S01CMH1d\\."]           = "201,94,255";
  cnA["S1C6H1L\\."]            = "178,255,229";
  cnA["S01C6H1L\\."]           = "178,255,229";
  cnA["S1C7H1L\\."]            = "63,135,171";
  cnA["S1C10H1L\\."]           = "94,94,255";
  cnA["S1C10H1-B\\."]          = "116,7,171";
  cnA["S1C10H1-C\\."]          = "255,10,173";
  cnA["S1C10H1-D\\."]          = "171,120,137";
  cnA["S1C10H2\\."]            = "171,99,63";
  cnA["S1C12H1L\\."]           = "62,171,7";
  cnA["S1C12H2\\."]            = "10,255,92";
  cnA["S1C16H1L\\."]           = "255,255,10";
  cnA["S1C1/5/19H1L\\."]       = "120,171,154";
  cnA["S1C5pH2\\."]            = "179,230,255";
  cnA["S1C5H2\\."]             = "179,230,255";
  cnA["S1C10/12H1d\\."]        = "255,204,179";
  cnA["S1CMH3d\\."]            = "255,204,179";
  cnA["S1C12H3d\\."]           = "171,7,116";
  cnA["S01C12H3d\\."]          = "171,7,116";
#  cnA["S1C12H4\\."]           = "171,171,63";
  cnA["S1C3H4-C\\."]           = "10,255,255";
  cnA["del_A_pos_70-GREYR.2"]="212,113,0";
  cnA["ref_A_pos_70-GREYR.2"]="14,138,30";
  cnA["del_27bp-GREYR.4"]="87,42,212";
  cnA["ref_27bp-GREYR.4"]="85,212,102";
  cnA["del_T_pos_150-YELLOW2.4"]="106,134,212";
  cnA["ref_T_pos_150-YELLOW2.4"]="83,138,119";
  cnA["del_7bp-YELLOW2.4"]="148,165,212";
  cnA["ref_7bp-YELLOW2.4"]="198,212,0";
  cnA["del_TTCT_pos_175-YELLOWG.4"]="47,14,138";
  cnA["ref_TTCT_pos_175-YELLOWG.4"]="212,42,155";
  cnA["ref_TC_pos_62-YELLOW2.4"]="55,138,66";
  cnA["del_A_pos_63-GREYB.5"]="69,87,138";
  cnA["ref_A_pos_63-GREYB.5"]="127,195,212";
  cnA["alt_pos_10_69-GREYB.5"]="165,148,212";
  cnA["ref_pos_10_69-GREYB.5"]="85,212,0";
  cnA["del_C_pos_159-GREYM.5"]="174,21,212";
  cnA["ref_C_pos_159-GREYM.5"]="212,42,65";
  cnA["del_14bp_pos_49_62-GREYZ.5"]="119,85,212";
  cnA["ref_14bp_pos_49_62-GREYZ.5"]="87,69,138";
  cnA["del_T_pos_135-YELLOW2.6"]="150,127,212";
  cnA["del_CCT_pos_60_62-YELLOW1.6"]="107,96,138";
  cnA["ref_CCT_pos_60_62-YELLOW1.6"]="0,212,141";
  cnA["del_AC_pos_60-GREYR.7"]="113,14,138";
  cnA["ref_AC_pos_60-GREYR.7"]="138,93,41";
  cnA["del_T_pos_135-YELLOW2.9"]="121,55,138";
  cnA["ref_T_pos_135-YELLOW2.9"]="134,106,212";
  cnA["del_A_pos_140-GREYP.10"]="195,127,212";
  cnA["ref_A_pos_140-GREYP.10"]="129,96,138";
  cnA["del_T_pos_135_and_pos_163-YELLOW2.10"]="0,110,138";
  cnA["ref_T_pos_135_and_pos_163-YELLOW2.10"]="138,14,96";
  cnA["del_T_pos_99-YELLOW3.11"]="123,212,64";
  cnA["ref_T_pos_99-YELLOW3.11"]="186,85,212";
  cnA["ins_SPL_W3J_GATC_pos_54-YELLOW2.10"]="212,106,176";
  cnA["ref_SPL_W3J_GATC_pos_54-YELLOW2.10"]="127,83,138";
  cnA["del_SPL_W4F_CCT_pos_61-YELLOW2.6"]="199,148,212";
  cnA["ref_SPL_W4F_CCT_pos_61-YELLOW2.6"]="0,56,212";
  cnA["alt_SPL_W4D_ATTTA_pos_59-YELLOW2.4"]="138,14,30";
  cnA["ref_SPL_W4D_CTTTTCC_pos_59-YELLOW2.4"]="64,212,162";

  # FEDOR updates:
  cnA["S3C1H2-F\\."]            = "0,255,0";
  cnA["S3C1H2-E\\."]            = "0,255,255";
  cnA["S2С2H2-C\\."]            = "0,255,0";
  cnA["S3C11H2-A\\."]           = "7,171,7";
  cnA["S3C11H2-B\\."]           = "0,255,255";
  cnA["S4C15H4\\."]             = "255,150,0";
  cnA["S5C5/19H8\\."]           = "0,0,255";
}

!/^#/ {
  chr_name ? targetName = chr_name : targetName = $1;
  queryName = $3;
# FEDOR: hmmer out is 1-based, bed is 0-based
  aliFrom   = $7-1;
  aliTo     = $8;
  envFrom   = $9;
  envTo     = $10;
  strand    = $12;
  score     = $14;
  description = $16;
  (description != "-" && description != "") ? description = " " $16 : description = "";

# skip SF monomers
# FEDOR: no. we need them
# if ((skipsf)&&(length(queryName) == 2)) { next }

# threshold score to length
  aliLength = abs(aliTo-aliFrom);
  if (+aliLength > 0) {
    if ((th)&&(sprintf("%.1f",score/aliLength)) < +th) { next }
  } else {
    next
  }

  if (coordsEnv) {
    if (strand == "+") {
      print targetName description, envFrom, envTo, queryName, sprintf("%.1f",score), strand, envFrom, envTo, getColor(queryName, cnA);
    } else {
      print targetName description, envTo, envFrom, queryName, sprintf("%.1f",score), strand, envTo, envFrom, getColor(queryName, cnA);
    }
  } else {
    if (strand == "+") {
      print targetName description, aliFrom, aliTo, queryName, sprintf("%.1f",score), strand, aliFrom, aliTo, getColor(queryName, cnA);
    } else {
      print targetName description, aliTo, aliFrom, queryName, sprintf("%.1f",score), strand, aliTo, aliFrom, getColor(queryName, cnA);
    }
  }
}

function abs(x){return ((x < 0.0) ? -x : x)}

function getColor(s, cnA,    color, i) {
  # here is assigned a color for each type
  color = "";
  for (i in cnA) {
    color = (match(s, i) ? cnA[i] : "85,85,85");
    if (color != "85,85,85") { break }
  }
  return color;
}
