# Stage 6C ASTRAL4/CASTLES-II provenance

Executable: `/Users/ytabatabaee/Desktop/ASTER/bin/astral4`
Executable SHA256: `2dc864161f7a80188f37fd6b39abc53c38aee548c20b534f274708b508a99b4c`
ASTER commit: `ddc3dc6b28b20cecfec18cede8c6247decc82b77`
Platform: `macOS-26.6.2-x86_64-i386-64bit`
Python: `3.12.11`
Threads: `2`
CASTLES-II status: the ASTRAL4 help output reports `NOW with integrated CASTLES-2`; `SULength` is the default substitution-per-site branch-length output.

## Help output

```text
Accurate Species TRee ALgorithm IV (ASTRAL-IV)
*** NOW with integrated CASTLES-2 ***
Version: v1.25.4.8
astral4 [ --ARG_NAME ARG_VALUE ... ] input_file
Abbr	Arg	Type	Default	Description
-C	--scoring	Preset		Scoring the full species tree file after `-c` without exploring other topologies (`-r 1 -s 0`)
-h	--help	Preset		Printing the help message
-R	--moreround	Preset		More rounds of placements and subsampling (`-r 16 -s 16`)
-o	--output	String	<standard output>	File name for the output species tree
-t	--thread	Integer	1	Number of threads

Advanced options:
-u	--support	Integer	1	Output support option (0: no branch or support, 1: length and support only, 2: detailed, 3: freqQuad.csv)
-w	--downweightrepeat	Float	1	The number of trees sampled for each locus
-l	--lambda	Float	0.5	Rate lambda of Yule process under which the species tree is modeled
-a	--mapping	String		A list of gene name to taxon name maps, each line contains one gene name followed by one taxon name separated by a space or tab
-v	--verbose	Integer	2	Level of logging (1: minimum, 2: normal)
-r	--round	Integer	4	Number of initial rounds of placements
-g	--guide	String		Newick file containing binary trees as guide trees
-c	--constraint	String		Newick file containing a binary species tree to place missing species on
-s	--subsample	Integer	4	Number of rounds of subsampling per exploration step
-i	--input	String	<last argument>	The input file
	--root	String		Root at the given species
	--seed	Integer	233	Seed for pseudorandomness
	--proportion	Float	0.25	Proportion of taxa in the subsample in naive algorithm
	--genelength	Float	1000	Average gene sequence length
	--sulengthtype	Integer	1	1: species tree branch length; 2: mean gene tree branch length
	--length	String	SULength	SULength: substitution-per-site unit; CULength: coalescent unit

Examples: 
astral4 -o output_file input_file
astral4 -i input_file -o output_file
```

## Commands

```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -C -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -c /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_fixed_background_topology.nwk -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/background_161.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_background.fixed_background_su.nwk
```
```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -C -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -c /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_fixed_background_topology.nwk -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/all_213.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_all.fixed_background_su.nwk
```
```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -C -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -c /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_fixed_background_topology.nwk -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/chr16_all_96.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_chr16_all.fixed_background_su.nwk
```
```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -C -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -c /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_fixed_background_topology.nwk -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/supergene_52.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_supergene.fixed_background_su.nwk
```
```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -R -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/background_161.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_background.free_su.nwk
```
```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -R -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/all_213.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_all.free_su.nwk
```
```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -R -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/chr16_all_96.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_chr16_all.free_su.nwk
```
```bash
/Users/ytabatabaee/Desktop/ASTER/bin/astral4 -R -u 2 -t 2 -a /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/processed/stage6c_astral4_mapping.txt -i /Users/ytabatabaee/Desktop/msrc-paper/data/fire_ants_chr16/intermediate/astral4_su/supergene_52.trees -o /Users/ytabatabaee/Desktop/msrc-paper/empirical/fire_ants_chr16/results/astral4_su/T_supergene.free_su.nwk
```
