#!/bin/bash -eux
set -o pipefail

set +x
eval $(/opt/homebrew/bin/brew shellenv)
set -x
PYENV_ROOT="${HOME}/.pyenv"
PATH="${PYENV_ROOT}/bin:${PATH}"
if command -v pyenv 1>/dev/null 2>&1; then
    eval "$(pyenv init --path)"
    eval "$(pyenv init -)"
fi
pyenv versions

script_dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
# the number of arguments is greater than one
if [ $# -ge 1 ]; then
  $script_dir/../../cmake-build-release/main_grouped_evo_ODE $@
else
  python $script_dir/get_args.py | xargs -t $script_dir/../../cmake-build-release/main_grouped_evo_ODE
fi
export PIPENV_PIPFILE=${script_dir}/../Pipfile
pipenv run python ${script_dir}/plot_timeseries.py
pipenv run python ${script_dir}/plot_histo.py
