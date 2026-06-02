#!/usr/bin/env bash
# Download and extract CellRanger 10.0.0
set -euo pipefail

INSTALL_DIR="${HOME}"
cd "${INSTALL_DIR}"

if [[ -x "${INSTALL_DIR}/cellranger-10.0.0/cellranger" ]]; then
  echo "CellRanger 10.0.0 already installed at ${INSTALL_DIR}/cellranger-10.0.0/"
  "${INSTALL_DIR}/cellranger-10.0.0/cellranger" --version
  exit 0
fi

echo "Downloading CellRanger 10.0.0..."
curl -o cellranger-10.0.0.tar.gz \
  "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1772424499&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=CePYH8zJN~BfQYjmq5A9K4rDEcsfusG8ak8Bax1-s7oM293F5XD582pSlDdE4v7l93D4DCZ~8-vQCfZiEU~ophZ7hDv7ddV21pu6YvyMBfVXexrSyFqRejaXbvqaQuYzUeiFTilmVSFbLISCeYS8I97waYV4-XGty5nXA~PPOsGsbnCiAcsgbyGOYmcAoo4nFIdRjVcgXxD0jWAgPG1R2oZ6OH6yuQ~~TyMZ0joBBJF7unjNnvjnT5Z41qtcAVpea0rSe1VS1WQ4zzTV47esSedaZmaQOORuxOs6ZtU~wKVr6qply05r3aOrVxT5ddowO0JrLRdSEVevjImhnU5iBQ__"

echo "Extracting..."
tar -xzf cellranger-10.0.0.tar.gz

echo "Cleaning up tarball..."
rm -f cellranger-10.0.0.tar.gz

echo "Verifying..."
"${INSTALL_DIR}/cellranger-10.0.0/cellranger" --version

echo "Done. Binary at: ${INSTALL_DIR}/cellranger-10.0.0/cellranger"
