# RNA3DGraphclust
This project aims to build a program to detect domains in RNA 3D structure derived from the PDB database, using several graph-based clustering algorithms.

### Prerequisites
PyMOL version 2.5 or later

### Installation
There are 2 ways to install the tool:

#### 1.  Using Docker:
```docker pull lequockhang/rna3dgraphclust ```

Then use this command below to run the docker image:
``` docker run -v `pwd`:/workdir/ lequockhang/rna3dgraphclust [options] ```

Whereas `` `pwd` `` is your full path to your working directory containing input file(s). Here is the full command for example:

``` docker run -v `pwd`:/workdir/ lequockhang/rna3dgraphclust -i Example.pdb -a G -v -p Output_pdb_name -j Output_json_name```

#### 2.  Using source code:
```git clone https://github.com/LeKhang97/RNA3DGraphclust```

From here, you can either build it globally or in a virtual environment:

##### 2.1 Build globally:
```pip3 install -r requirements.txt```

##### 2.2 Build in a virtual environment:
```make```

Then you can execute the program in virtual environment by:
```./venv/bin/python3 RNA3DGraphclust.py [options]```

### Usage
You can execute the program by:<br/>
```python3 RNA3DGraphclust.py -i infile -v -a M -o path_for_output  ```

Type ```python3 RNA3DGraphclust.py -h``` for more information of the usage:
```
positional arguments:
  {G,M,L,C,H}
    G                   Arguments for Girvan-Newman clustering algorithm
    M                   Arguments for Markov clustering algorithm
    L                   Arguments for Louvain clustering algorithm
    C                   Arguments for Clauset-Newman-Moore clustering algorithm
    H                   Arguments for Hierachical-based clustering algorithm

options:
  -h, --help            show this help message and exit
  -v, --verbose         verbose mode.
  -i INPUT, --input INPUT
                        input file. Must be in pdb format.
  -at, --atom_type      Atom types to be considered in the analysis. The default is C3.

  -t THRESHOLD, --threshold THRESHOLD
                        Lower threshold for sequence length.
  -o OUTPATH, --outpath OUTPATH
                        path of output for json and pdb files. If not specified, the output will be saved in the current directory.

  -p PDB, --pdb PDB
                        output filename in PDB format.

  -j JSON, --json JSON
                        output filename in JSON format.

  -a {G,M,L,C,H}, --algorithm {G,M,L,C,H}
                        Clustering algorithm. Either: G (Girvan-Newman, default); M (Markov); L (Louvain); C (Clauset-Newman-Moore); H (Hierachical-based))
```

- Each algorithm has its default parameters. For example, if you want to check the Markov clustering, type ```python3 ./RNA3DGraphclust.py M -h ``` for details. You can also change the parameters, in this case is the epsilon (-e), by following: <br>
``` python3 RNA3DGraphclust.py -i infile -v -a M -o . M -e 5```

### Notes
- All parameters specified for each algorithm can only be accessed after typing its abbreviation (besides the option -a Algorithm);
- The input must be in pdb format;
- There are 2 output files if the output option is chosen. One file is the **JSON file**, which contains the coordinate, the residue number and the label of clusters. The other file contains the **command line for PyMOL GUI** to generate the clusters, which has the same name as the JSON file with the suffix '_pymolcmd.pml' (name ```outfile_pymolcmd.pml``` for example). You can either run it from terminal by this command:<br>
`pymol outfile_pymolcmd.pml`
<br/> if you already created alias for PyMOL. Or from PyMOL command line, you can try: <br/>
```@outfile_pymolcmd ```
