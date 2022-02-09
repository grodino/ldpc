# Error correcting codes toolbox

This repository contains functions to generate, decode and simulate LDPC (Low Density Sparse Graph) codes for research purposes. It has been created as a mini-project for the 2022 [Master SAR](https://www.universite-paris-saclay.fr/formation/master/electronique-energie-electrique-automatique/m2-advanced-wireless-communications-systems) *Coding theory* course with [Prof. Antoine Berthet](https://l2s.centralesupelec.fr/u/berthet-antoine/).
You can do whatever with this code as long as you credit me.

## Installation

There are two options to run the code: either use a precompiled binary or build it with [Cargo](https://www.rust-lang.org/tools/install)

### Use precompiled binary
Go to the [release page](https://github.com/grodino/ldpc/releases) of the repository and download the binary 
corresponding to your system.

#### Windows
Just run the binary : 
```shell script
ldpc.exe --help
```

#### Linux
Set the file permissions as executable : `chmod +x ldpc-linux` and run it :
```shell script
./ldpc-linux --help
```

#### MacOs
A binary is provided but not tested

### Install via Cargo 
This is the most reliable method. Install [Rust](https://www.rust-lang.org/tools/install) on your machine, clone the 
repository and run cargo.

Once Rust is installed, run :
```shell script
git clone https://github.com/grodino/ldpc.git
cd ldpc
cargo run -- --help
```

> The `--` after `cargo run` are here to pass the command line arguments to `ldpc` after it has been built by `cargo`.


## Report
Here are the commands used to create the figures in the report:

### Gallager codes
```bash
ldpc --snr="-3,7,20" -d 3 gallager > scripts/results/gallager.json
python scripts/experiments.py scripts/results/gallager.json --save
```

### Regular Mac Kay & Neal codes
```bash
ldpc --snr="-3,7,20" -d 2 mac-kay -r > scripts/results/mac-kay_regular.json
python scripts/experiments.py scripts/results/mac-kay_regular.json --save
```

### Irregular Mac Kay & Neal codes
```bash
ldpc --snr="-3,7,20" -d 2 mac-kay > scripts/results/mac-kay.json
python scripts/experiments.py scripts/results/mac-kay.json --save
```

### Irregular Mac Kay & Neal codes
```bash
ldpc --snr="-3,7,20" -d 2 mac-kay -f > scripts/results/mac-kay_full-rank.json
python scripts/experiments.py scripts/results/mac-kay_full-rank.json --save
```