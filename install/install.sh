#!/usr/bin/env bash


SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
PATH_SCRIPT="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

PATH_BASE="${PATH_SCRIPT}/.."
PATH_BIN="${PATH_BASE}/bin"


PATH_BASH_DRIVER=${PATH_BASE}/code/bash/driver
PATH_PYTHON_DRIVER=${PATH_BASE}/code/python/driver
PATH_PERL_DRIVER=${PATH_BASE}/code/perl/driver


function create_exec_script_with_cmd() {

  local cmd="$1"
  local script_name="$2"

  echo "#!/usr/bin/env bash" > ${PATH_BIN}/${script_name}
  # echo "DIR=\"\$( cd \"\$( dirname \"\${BASH_SOURCE[0]}\" )\" >/dev/null && pwd )\"" >> ./bin/${script_name}
  echo "$cmd" >> ${PATH_BIN}/${script_name}

  chmod +x ${PATH_BIN}/${script_name}
}

function create_name_for_exec_script() {
  local name_of_script_to_execute="$1"

  local script_name_no_ext=$(get_filename_without_path_or_extension "${name_of_script_to_execute}")
  local script_ext=$(get_file_extension "${name_of_script_to_execute}")

  local name_of_exec="${script_name_no_ext}_${script_ext}.sh"

  echo "$name_of_exec"
}

function get_interpreter_from_extension() {
  case ${script_ext} in
    "sh")
      env_type="source"
      ;;
    "py")
      env_type="python"
      ;;
    "pl")
      env_type="perl"
      ;;
    *)
      env_type=
      ;;
    esac

    echo "${env_type}"

}

function create_exec_script_for_directory() {

  local driver_dir="$1";

  # for every script in bash driver, create a separate bin script for it
  for script in ${driver_dir}/*; do

    [ -f "$script" ] || continue                                # make sure it's a file
    [[ ! "${script}" =~ .*-template\.* ]] || continue           # ignore templates

    local name_of_exec=$(create_name_for_exec_script ${script})
    local script_ext=$(get_file_extension "${script}")

    interpreter=$(get_interpreter_from_extension ${script_ext})

    if [[ -z ${interpreter} ]]; then
      >&2 echo "Script has unknown extension. Ignoring... $script"
      continue
    fi



    local exec_cmd="$interpreter ${ABS_PATH_SCRIPT}/$script \"\$@\""

    echo "Creating: $name_of_exec"

    
    create_exec_script_with_cmd "${exec_cmd}" "${name_of_exec}"

  done

}

function get_filename_without_path_or_extension() {
  local fullfilename="$1"
  local filename=$(basename "$fullfilename")
  local fname="${filename%.*}"
  echo "$fname"
}

function get_file_extension() {
  local fullfilename="$1"
  local filename=$(basename "$fullfilename")
  local ext="${filename##*.}"
  echo "$ext"
}

#create_exec_script_for_directory $PATH_BASH_DRIVER
create_exec_script_for_directory $PATH_PYTHON_DRIVER
#create_exec_script_for_directory $PATH_PERL_DRIVER
