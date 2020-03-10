#!/usr/bin/env bash
rm *.zip
zip PRSice_linux.zip PRSice.R TOY* PRSice_linux
zip PRSice_mac.zip PRSice.R TOY* PRSice_mac
zip PRSice_win64.zip PRSice.R TOY* PRSice_win64.exe
zip PRSice_win32.zip PRSice.R TOY* PRSice_win32.exe
