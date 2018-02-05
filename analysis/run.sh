#!/usr/bin/env bash

if [[ "$(uname)" == "Darwin" ]];
then
  script=$(greadlink -f $0)
else
  script=$(readlink -f $0)
fi
script_path=`dirname $script`

"${script_path}/1_normalize_data.R"
"${script_path}/1_normalize_data_plot.R"

"${script_path}/2_model_selection.R"

"${script_path}/3_pan_viral_hit_selection.R"
"${script_path}/3_pan_viral_hit_selection_plot.R"

"${script_path}/4_stability_selection.R"
"${script_path}/4_stability_selection_plot.R"

"${script_path}/5_network_analysis.R"
"${script_path}/5_network_analysis_table.R"
