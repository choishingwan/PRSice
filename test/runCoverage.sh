#!/bin/sh

# ${1} is the directory containing the .gcno files (%{buildDir} in Qt Creator)
# ${2} is the unit test executable
# ${3} is the directory for the unit test output
LCOV=lcov
GENHTML=genhtml
#BROWSER=Firefox 

SRC_DIR="${1}"
UNIT_TEST="${2}"
TEST_DIR="${3}"

# Generate our base line info
"${LCOV}" -d "${SRC_DIR}" -z
"${LCOV}" --rc lcov_branch_coverage=1 --capture -i -b "${SRC_DIR}/../"  -d "${SRC_DIR}" -o "${SRC_DIR}/baseline.info" > log

"${UNIT_TEST}" "${3}"

HTML_RESULTS="${1}/html"

mkdir -p ${HTML_RESULTS}

# generate our coverage info
"${LCOV}" --rc lcov_branch_coverage=1  --capture -b "${SRC_DIR}/../" -d "${SRC_DIR}" -o "${SRC_DIR}/coverage.info" > log

# Combine the two info
"${LCOV}" -a "${SRC_DIR}/baseline.info" -a "${SRC_DIR}/coverage.info" -o "${SRC_DIR}/coverage-combined.info" >log

# remove some paths
"${LCOV}" -r "${SRC_DIR}/coverage-combined.info" "*gtest*" "*lib*" "*Qt*.framework*" "*Xcode.app*" "*.moc" "*moc_*.cpp" "*/test/*" -o "${SRC_DIR}/coverage-filtered.info"  >log

#"${LCOV}" -r "${SRC_DIR}/coverage.info" "*gtest*" 

# generate our HTML
"${GENHTML}" --branch-coverage -o "${HTML_RESULTS}" "${SRC_DIR}/coverage-filtered.info" >log

# reset our counts
"${LCOV}" -d "${SRC_DIR}" -z

# open in browser and bring to front
#open -a "${BROWSER}" "${HTML_RESULTS}/index.html"
