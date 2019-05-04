#!/bin/sh

# ${1} is the directory containing the .gcno files (%{buildDir} in Qt Creator)
# ${2} is the unit test executable
# ${3} is the directory for the unit test output
LCOV=lcov
GENHTML=genhtml
BROWSER=Firefox 

SRC_DIR="${1}"
UNIT_TEST="${2}"
TEST_DIR="${3}"

# Generate our base line info

"${LCOV}" -z -c -i -b "${SRC_DIR}/../"  -d "${SRC_DIR}" -o "${SRC_DIR}/baseline.info"

"${UNIT_TEST}" "${3}"

HTML_RESULTS="${1}/html"

mkdir -p ${HTML_RESULTS}

# generate our coverage info
"${LCOV}" -b "${SRC_DIR}/../" -d "${SRC_DIR}" -c -o "${SRC_DIR}/coverage.info"

# Combine the two info
"${LCOV}" -a "${SRC_DIR}/baseline.info" -a "${SRC_DIR}/coverage.info" -o "${SRC_DIR}/coverage-combined.info"

# remove some paths
"${LCOV}" -r "${SRC_DIR}/coverage-combined.info" "*gtest*" "*lib*" "*Qt*.framework*" "*Xcode.app*" "*.moc" "*moc_*.cpp" "*/test/*" -o "${SRC_DIR}/coverage-filtered.info"

#"${LCOV}" -r "${SRC_DIR}/coverage.info" "*gtest*"

# generate our HTML
"${GENHTML}" -o "${HTML_RESULTS}" "${SRC_DIR}/coverage-filtered.info"

# reset our counts
"${LCOV}" -d "${SRC_DIR}" -z

# open in browser and bring to front
open -a "${BROWSER}" "${HTML_RESULTS}/index.html"
