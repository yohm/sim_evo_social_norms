#!/bin/bash -eux

script_dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
# the number of arguments is greater than one
if [ $# -ge 1 ]; then
  $script_dir/../cmake-build-release/main_grouped_evo $@
else
  python $script_dir/get_args.py | xargs -t $script_dir/../cmake-build-release/main_grouped_evo
fi
cat histo_norms.dat
export PIPENV_PIPFILE=${script_dir}/Pipfile
pipenv run python ${script_dir}/plot_timeseries.py
pipenv run python ${script_dir}/plot_histo_norms.py
