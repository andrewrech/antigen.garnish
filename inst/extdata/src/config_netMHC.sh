#!/bin/bash

# from the directory containing the four netMHC tarballs and the antigen.garnish
# data directory (downloaded from the one-line installation script)
# copy the tarballs into the antigen.garnish directory
NET_MHC_DIR=/
ANTIGEN_GARNISH_DIR=/root/antigen.garnish

cd "$NET_MHC_DIR" || return 1

mkdir -p "$ANTIGEN_GARNISH_DIR/netMHC" || return 1

find . -name "netMHC*.tar.gz" -exec tar xvzf {} -C "$ANTIGEN_GARNISH_DIR/netMHC" \;

chown "$USER" "$ANTIGEN_GARNISH_DIR/netMHC"
chmod 700 -R "$ANTIGEN_GARNISH_DIR/netMHC"

# antigen.garnish is ready to run in R!
