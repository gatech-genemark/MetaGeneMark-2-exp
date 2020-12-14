
# get path of config script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


export base=$DIR
export bin=$base/bin
export bin_external=$base/bin_external

export data=$base/data
export metadata=$base/metadata

export lists=$base/lists

export tmp=$base/tmp               # it is suggested that runs be performed in tmp directory (under some subdirectory)
export config=$base/config
export runs=$base/runs
export code=$base/code

export libpython=$code/python/lib
export driverpython=$code/python/driver

