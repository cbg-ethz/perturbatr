#!/usr/bin/env bash

if [[ "$(uname)" == "Darwin" ]];
then
  script=$(greadlink -f $0)
else
  script=$(readlink -f $0)
fi
script_path=`dirname $script`

"${script_path}/1a_normalize_data.R"
"${script_path}/1b_normalize_data_plot.R"

"${script_path}/2_model_selection.R"

"${script_path}/3a_pan_viral_hit_selection.R"
"${script_path}/3b_pan_viral_hit_selection_plot.R"

"${script_path}/4a_stability_selection.R"
"${script_path}/4b_stability_selection_plot.R"

"${script_path}/4c_predictive_model.R"
"${script_path}/4d_predictive_model_plot.R"

"${script_path}/5a_network_analysis.R"
"${script_path}/5b_network_analysis_table.R"
