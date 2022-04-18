#!/bin/bash
# This is a modified version of the pre-commit git script that checks for any mismatched
# md/Rmd readmes (including in subdirectories)

README_RMD=($(git diff --cached --name-only | grep -Ei 'README\.Rmd$'))
README_MD=($(git diff --cached --name-only | grep -Ei 'README\.md$'))

MSG="use 'git commit --no-verify' to override this check"

# No readme files
if [[ ${#README_RMD[@]} == 0 && ${#README_MD[@]} == 0 ]]; then
  exit 0
fi

# No readme files
if [[ ${#README_RMD[@]} -lt ${#README_MD[@]} ]]; then
  echo -e "There are more README.md files committed than README.Rmd files.\nPlease stage one Rmd for each md\n$MSG"
  exit 1
fi


# Checks to see if the RMD files and MD files are matched
function compare_rmd {
  local RMD=$1
  local DIR=$(dirname $RMD)
  
  if [[ $DIR/README.Rmd -nt $DIR/README.md ]]; then
    echo -e "$DIR/README.md is out of date; please re-knit $DIR/README.Rmd\n$MSG"
    exit 1
  fi
  
  if [[ $DIR == "." ]]; then 
    # Special case for base readme, since the general case fails
    if [[ "${README_MD[@]}" != *" README"* && "${README_MD[@]}" != "README"* ]]; then 
    # Account for spaces/starting position
      echo -e "README.Rmd and README.md should be both staged\n$MSG"
      exit 1
    fi
  elif [[ "${README_MD[@]}"  != *"${DIR}/README"* ]]; then # Both are present
    
      # echo "${README_MD[@]}"
      # echo "Dir: $DIR"
      # echo "RMD: $RMD"
      echo -e "$DIR/README.Rmd and $DIR/README.md should be both staged\n$MSG"
      exit 1
  fi
}

for rmd in "${README_RMD[@]}"; do
  compare_rmd $rmd
done