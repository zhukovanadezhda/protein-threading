# Protein Structure Prediction Using Threading


This project focuses on protein structure prediction using the [threading (fold recognition) method](https://en.wikipedia.org/wiki/Threading_(protein_sequence)). Threading is a computational approach that aligns a protein sequence with known template structures to predict its 3D structure, even when sequence similarity is low. This method identifies candidate templates by considering structural similarities, such as predicted secondary structures and solvent accessibility. By mapping the sequence onto structurally similar templates, threading can accurately predict protein folds, making it a valuable tool in bioinformatics, medicine, and biotechnology. The algorithm is inspired by David Jones's THREADER *(Computational Methods in Molecular Biology, Chapter 13, 1998)*.

<p align="center">
  <img align="center" width="800px" 
    src="doc/assets/fold-recognition.png"
    alt="protein">
</p>
<p align="center">
  <i>
    Fig. 1: Protein fold recognition summary.
  </i>
</p>


## 🔄Installation

### Clone the repository

```bash
git clone git@github.com:zhukovanadezhda/protein-threading.git
cd protein-threading
```
### Setup the conda environment

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html). Create the `protein-threading` conda environment:

```bash
conda env create -f environment.yml
```

### Load the environment

```bash
conda activate protein-threading
```

> 💡**Note:** To deactivate an active environment, use:
> ```bash
> conda deactivate
> ```

## 🧑‍💻️ Usage

> 💡**Note:** Before running the program, ensure the `src/config.py` file is properly configured to set up your working directories.  
> - **To run an example:** Select one of the provided example directories (e.g., `data/example1`).  
> - **To use custom data:** Create your own directory (e.g., `data/your_dir`) with subdirectories for `sequences` and `structures`. Place your data in these folders and update the paths in `src/config.py`.
>   
> Example `src/config.py` modification:
> ```python
> # Directory paths
> TEMPLATES_DIR = 'data/your_dir/structures/'
> SEQUENCES_DIR = 'data/your_dir/sequences/'
> ```

To run the program, use the following command:

```python
python src/main.py [-h] [--sequences SEQUENCES] [--templates TEMPLATES] [--output_file OUTPUT_FILE] \
                   [--jobs JOBS] [--dry_run] [--verbose]
```

## ⚙️Arguments

| Argument                  | Description                                                   | Default           |
|:-------------------------:|---------------------------------------------------------------|-------------------|
| `-h`                      | Show a help message and exit.                                 |                   |
| `--sequences`             | Comma-separated list of sequence filenames (`.fasta` format). | All files from `SEQUENCES_DIR` from `src/config.py`. |
| `--templates`             | Comma-separated list of template filenames (`.pdb` format).   | All files from `TEMPLATES_DIR` from `src/config.py`. |
| `--gap_score`             | The gap penalty.                                              | `0`|
| `--output_file`           | Name of the output CSV file.                                  | `results/energy_scores.csv`|
| `--jobs`                  | Number of parallel jobs to run.                               | All cores         |
| `--dry_run`               | If set, only log actions without processing.                  | `False` (not set)   |
| `--verbose`               | If set, verbose output is enabled.                            | `False` (not set)   |
| `--print_alignments`      | If set, the alignments are printed.                           | `False` (not set)   |

<p align="center">
  <i>
    Table 1: Program parameters.
  </i>
</p>

## 🎁Examples

### ☝️Example 1: Small proteins (<50 amino acids)  
**Estimated execution time: Less than 35 minutes**

#### Input

> `src/config.py` :
> ```python
> # Directory paths
> TEMPLATES_DIR = 'data/example1/structures/'
> SEQUENCES_DIR = 'data/example1/sequences/'
> ```

```python
python src/main.py --sequences 1CRN.fasta,1L2Y.fasta,1VII.fasta,5AWL.fasta --gap_score 0.2 \
                   --output_file results/example1_result.csv --print_alignments
```

#### Output
```bash
2024-09-17 08:10:16,722 - INFO - Loading DOPE score data...
2024-09-17 08:11:47,306 - INFO - Processing sequences and templates...
2024-09-17 08:11:47,306 - INFO - Processing sequence 5AWL.fasta, length: 10
2024-09-17 08:11:50,621 - INFO - Processing template 5awl.pdb with 10 residues.
2024-09-17 08:11:50,444 - INFO - Processing template 1le0.pdb with 12 residues.
  0  1  2  3  4  5  6  7  8  9
  |  |  |  |  |  |  |  |  |  |
  Y  Y  Y  D  P  E  T  G  T  W
2024-09-17 08:11:57,940 - INFO - Processed template 5awl.pdb. Energy score: -32.68
  0  1  2  3  4  5  6  7  8  9 10 11
  |  |  |  |  |  |  |  |     |  |   
  Y  Y  Y  D  P  E  T  G  -  T  W  -
2024-09-17 08:12:00,712 - INFO - Processed template 1le0.pdb. Energy score: -28.12
...
```


#### Results

<div style="text-align:center;">
  <table align="center">
    <thead>
      <tr>
        <th style="text-align:center;"> </th>
        <th style="text-align:center; width: 150px;">1crn.pdb <br> <img src="doc/assets/1crn.png" align="center" alt="1crn" width="70"></th>
        <th style="text-align:center; width: 150px;">1l2y.pdb <br> <img src="doc/assets/1l2y.png" align="center" alt="1l2y" width="70"></th>
        <th style="text-align:center; width: 150px;">1le0.pdb <br> <img src="doc/assets/1le0.png" align="center" alt="1le0" width="70"></th>
        <th style="text-align:center; width: 150px;">1le1.pdb <br> <img src="doc/assets/1le1.png" align="center" alt="1le1" width="70"></th>
        <th style="text-align:center; width: 150px;">1vii.pdb <br> <img src="doc/assets/1vii.png" align="center" alt="1vii" width="70"></th>
        <th style="text-align:center; width: 150px;">5awl.pdb <br> <img src="doc/assets/5awl.png" align="center" alt="5awl" width="70"></th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td> </td>
        <td><p align="center">46 aa</p></td>
        <td><p align="center">20 aa</p></td>
        <td><p align="center">12 aa</p></td>
        <td><p align="center">12 aa</p></td>
        <td><p align="center">36 aa</p></td>
        <td><p align="center">10 aa</p></td>
      </tr>
      <tr>
        <td><p align="center">1CRN.fasta</p><img src="doc/assets/1crn.png" alt="1CRN" width="70"></td>
        <td><p align="center"><i><b>-283.99</b></i></p></td>
        <td><p align="center">-35.37</p></td>
        <td><p align="center">-10.5</p></td>
        <td><p align="center">-13.05</p></td>
        <td><p align="center">-162.04</p></td>
        <td><p align="center">-9.02</p></td>
      </tr>
      <tr>
        <td><p align="center">1L2Y.fasta</p><img src="doc/assets/1l2y.png" alt="1L2Y" width="70"></td>
        <td><p align="center">-14.77</p></td>
        <td><p align="center"><b><i>-67.4</p></b></i></td>
        <td><p align="center">-19.18</p></td>
        <td><p align="center">-20.23</p></td>
        <td><p align="center">-41.71</p></td>
        <td><p align="center">-17.53</p></td>
      </tr>
      <tr>
        <td><p align="center">1VII.fasta</p><img src="doc/assets/1vii.png" alt="1VII" width="70"></td>
        <td><p align="center">-157.37</p></td>
        <td><p align="center">-24.16</p></td>
        <td><p align="center">7.61</p></td>
        <td><p align="center">7.23</p></td>
        <td><p align="center"><b><i>-169.07</b></i></p></td>
        <td><p align="center">6.95</p></td>
      </tr>
      <tr>
        <td><p align="center">5AWL.fasta</p><img src="doc/assets/5awl.png" alt="5AWL" width="70"></td>
        <td><p align="center">6.73</p></td>
        <td><p align="center">-19.18</p></td>
        <td><p align="center">-28.12</p></td>
        <td><p align="center">-28.32</p></td>
        <td><p align="center">3.62</p></td>
        <td><p align="center"><b><i>-32.68</i></b></p></td>
      </tr>
    </tbody>
  </table>
</div>

<p align="center">
  <i>
    Table 2: Summary of results from the first example.<br>
    "aa" stands for amino acids; bold indicates the best score, italics indicate the correct match.
  </i>
</p>


### ✌️Example 2: Medium-sized proteins (70-120 amino acids)  
**Estimated execution time: Over 16 hours**

#### Input

> `src/config.py` :
> ```python
> # Directory paths
> TEMPLATES_DIR = 'data/example2/structures/'
> SEQUENCES_DIR = 'data/example2/sequences/'
> ```


```python
python src/main.py --sequences 1E68.fasta,3E8V.fasta --gap_score 0.1 --output_file results/example2_result.csv
```

#### Results

<div style="text-align:center;">
  <table align="center">
    <thead>
      <tr>
        <th style="text-align:center;"> </th>
        <th style="text-align:center; width: 150px;">1e68.pdb <br> <img src="doc/assets/1e68.png" align="center" alt="1e68" width="70"></th>
        <th style="text-align:center; width: 150px;">1ubq.pdb <br> <img src="doc/assets/1ubq.png" align="center" alt="1ubq" width="70"></th>
        <th style="text-align:center; width: 150px;">3zow.pdb <br> <img src="doc/assets/3zow.png" align="center" alt="3zow" width="70"></th>
        <th style="text-align:center; width: 150px;">3e8v.pdb <br> <img src="doc/assets/3e8v.png" align="center" alt="3e8v" width="70"></th>
        <th style="text-align:center; width: 150px;">1tit.pdb <br> <img src="doc/assets/1tit.png" align="center" alt="1tit" width="70"></th>
        <th style="text-align:center; width: 150px;">1tvd.pdb <br> <img src="doc/assets/1tvd.png" align="center" alt="1tvd" width="70"></th>
        <th style="text-align:center; width: 150px;">3zbv.pdb <br> <img src="doc/assets/3zbv.png" align="center" alt="3zbv" width="70"></th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td><p align="center"></p></td>
        <td><p align="center">70 aa</p></td>
        <td><p align="center">76 aa</p></td>
        <td><p align="center">81 aa</p></td>
        <td><p align="center">82 aa</p></td>
        <td><p align="center">89 aa</p></td>
        <td><p align="center">116 aa</p></td>
        <td><p align="center">118 aa</p></td>
      </tr>
      <tr>
        <td><p align="center">1E68.fasta <br> <img src="doc/assets/1e68.png" align="center" alt="1e68" width="70"></p></td>
        <td><p align="center"><i><b>-490.53</b></i></p></td>
        <td><p align="center">-466.0</p></td>
        <td><p align="center">-472.0</p></td>
        <td><p align="center">-480.22</p></td>
        <td><p align="center">-412.13</p></td>
        <td><p align="center">-316.45</p></td>
        <td><p align="center">-312.76</p></td>
      </tr>
      <td><p align="center">3E8V.fasta <br> <img src="doc/assets/3e8v.png" align="center" alt="3e8v" width="70"></p></td>
        <td><p align="center">-426.15</p></td>
        <td><p align="center">-479.47</p></td>
        <td><p align="center">-527.57</p></td>
        <td><p align="center"><i><b>-558.23</b></i></p></td>
        <td><p align="center">-493.82</p></td>
        <td><p align="center">-415.08</p></td>
        <td><p align="center">-407.47</p></td>
      </tr>
    </tbody>
  </table>
</div>
          
<p align="center">
  <i>
    Table 3: Summary of results from the second example.<br>
    "aa" stands for amino acids; bold indicates the best score, italics indicate the correct match.
  </i>
</p>

## 📊 Evaluate results significance

To determine whether the results obtained are significantly different from those expected by chance, a statistical approach based on the calculation of z-scores is used. The energy scores observed for each alignment can be compared to a distribution generated by a random shuffle of the sequences. The z-score quantifies the deviation from the mean of the random scores: a z-score lower than `-1.96` indicates, with 95\% confidence, a significant structural match, while a z-score higher than `1.96` signals an alignment significantly worse than chance. A z-score close to zero suggests the absence of a significant deviation.

### Usage:

> 💡**Note:** Before running the program, ensure the `src/config.py` file is properly configured with the correct paths for your working directories.  
> The working directories should correspond to those where the sequence and template files referenced in the energy scores file are stored.


```python
python src/evaluate_significance.py --input_csv <path_to_input_csv> --output_file <path_to_output_csv> \
                                    --gap_score <gap_score> --n_shuffles <number_of_shuffles>
```

### Requirements:
`input_csv` - the input csv file with calculated energy scores, the file should be in the format produced by the main script described above.

### Example:

#### Input 
> `src/config.py` :
> ```python
> # Directory paths
> TEMPLATES_DIR = 'data/example1/structures/'
> SEQUENCES_DIR = 'data/example1/sequences/'
> ```

```python
python src/evaluate_significance.py --input_csv results/example1_result.csv \
                                    --output_file results/shuffle_example1.csv \
                                    --n_shuffle 10 --gap_score 0.2
```

#### Output

```bash
WARNING - The number of shuffles is less than 30. Shapiro-Wilk test will be performed to check for the normality of the shuffled scores.
WARNING - The distribution of shuffled energy scores for sequence 5AWL.fasta is not normally distributed (p-value = 0.0000).
...
```

#### Result

<div style="text-align:center;">
  <table align="center">
    <tr>
      <th></th>
      <th><p align="center"><br>1crn.pdb</br></p></th>
      <th><p align="center"><br>1l2y.pdb</br></p></th>
      <th><p align="center"><br>1le0.pdb</br></p></th>
      <th><p align="center"><br>1le1.pdb</br></p></th>
      <th><p align="center"><br>1vii.pdb</br></p></th>
      <th><p align="center"><br>5awl.pdb</br></p></th>
    </tr>
    <tr>
      <td><p align="center"><br>1CRN.fasta</br></p></td>
      <td><p align="center"><br>-2.13</br></p></td>
      <td><p align="center">0.68</p></td>
      <td><p align="center">0.28</p></td>
      <td><p align="center">-0.91</p></td>
      <td><p align="center">1.22</p></td>
      <td><p align="center">-0.39</p></td>
    </tr>
    <tr>
      <td><p align="center"><br>1L2Y.fasta</br></p></td>
      <td><p align="center">1.12</p></td>
      <td><p align="center"><br>-2.49</br></p></td>
      <td><p align="center">2.96</p></td>
      <td><p align="center">2.62</p></td>
      <td><p align="center">1.69</p></td>
      <td><p align="center">1.66</p></td>
    </tr>
    <tr>
      <td><p align="center"><br>1VII.fasta</br></td>
      <td><p align="center">-0.16</p></td>
      <td><p align="center">0.18</p></td>
      <td><p align="center">-1.05</p></td>
      <td><p align="center">-1.40</p></td>
      <td><p align="center">-1.48</p></td>
      <td><p align="center"><br>-2.71</br></p></td>
    </tr>
    <tr>
      <td><p align="center"><br>5AWL.fasta</br></td>
      <td><p align="center">-0.97</p></td>
      <td><p align="center">-0.16</p></td>
      <td><p align="center">-1.32</p></td>
      <td><p align="center">-1.53</p></td>
      <td><p align="center">-1.19</p></td>
      <td><p align="center">-1.56</p></td>
    </tr>
  </table>
</div>
<p align="center">
  <i>
    Table 4: Summary of z-scores for alignments between different structures in example 1.
  </i>
</p>

> ⚠️ **Warning:** The z-scores may not be interpretable as the shuffled energy scores are not normally distributed. This is merely an example. To perform a meaningful significance evaluation, increase the number of shuffles.

## 🔧 Test Different Gap Penalties 

To explore the impact of gap penalties on sequence-structure matching performance, you can use the `test_gaps.py` script. This script evaluates several gap penalty values to determine the optimal setting for the algorithm. Performance is assessed based on the following scoring system:

- **Correct Structure Prediction**: If the algorithm accurately predicts the structure of the sequence, it earns **2 points**.
- **Similar Structure Bonus**: If a similar structure appears among the top two predictions, the algorithm earns an additional **1 point**.
- Thus, the maximum score for each sequence is **3 points**, and the final score is normalized by dividing it by the maximum possible score.

### Usage Instructions

> 💡 **Note**: Before running the program, ensure that `src/config.py` is properly configured. Set the paths for your working directories, the gap penalties to test, and similar structures to evaluate performance.
>> Example modification to `src/config.py`:
>>
>> ```python
>> # Gap penalties to test 
>> gap_scores = [0, 0.1, 0.2] 
>>
>> # Proteins with similar structures 
>> homolog_pairs = { 
>>     'A.fasta': ['b.pdb'], 
>>     'B.fasta': ['a.pdb', 'c.pdb'], 
>>     'C.fasta': ['b.pdb'] 
>> } 
>> ``` 


### Usage:

```python
python src/test_gaps.py [--program_path PROGRAM_PATH] [--output_dir OUTPUT_DIR] [--result_file RESULT_FILE]
```

### Example:

#### Input 
> `src/config.py` :
> ```python
> # Gap scores to test
> gap_scores = [0, 0.1, 0.2, 0.3, 0.5, 1, 2, 5]
>
> # Proteins with similar structures
> homolog_pairs = {
>     '5AWL.fasta': ['1l2y.pdb', '1vii.pdb', '1crn.pdb'],
>     '1VII.fasta': ['1l2y.pdb', '1crn.pdb'],
>     '1L2Y.fasta': ['1vii.pdb', '1crn.pdb'],
>     '1CRN.fasta': ['1l2y.pdb','1vii.pdb', '1crn.pdb', '1le0.pdb', '1le1.pdb']
> }
> ```

```python
python src/test_gaps.py --output_dir results/gaps_test --result_file results/performance.csv
```

#### Result

<div style="text-align:center;">
  <table align="center">
    <tr>
      <th><p align="center"><br>Gap Score</br></p></th>
      <th><p align="center"><br>Performance</br></p></th>
      <th><p align="center"><br>Correctly Guessed</br></p></th>
      <th><p align="center"><br>Similar Structure</br></p></th>
    </tr>
    <tr>
      <td><p align="center">0.0</p></td>
      <td><p align="center">0.5</p></td>
      <td><p align="center">1</p></td>
      <td><p align="center">4</p></td>
    </tr>
    <tr>
      <td><p align="center">0.1</p></td>
      <td><p align="center">0.58</p></td>
      <td><p align="center">2</p></td>
      <td><p align="center">3</p></td>
    </tr>
    <tr>
      <td><p align="center">0.2</p></td>
      <td><p align="center">0.92</p></td>
      <td><p align="center">4</p></td>
      <td><p align="center">3</p></td>
    </tr>
    <tr>
      <td><p align="center">0.3</p></td>
      <td><p align="center">0.92</p></td>
      <td><p align="center">4</p></td>
      <td><p align="center">3</p></td>
    </tr>
    <tr>
      <td><p align="center">0.5</p></td>
      <td><p align="center">0.83</p></td>
      <td><p align="center">4</p></td>
      <td><p align="center">2</p></td>
    </tr>
    <tr>
      <td><p align="center">1.0</p></td>
      <td><p align="center">0.75</p></td>
      <td><p align="center">4</p></td>
      <td><p align="center">1</p></td>
    </tr>
    <tr>
      <td><p align="center">2.0</p></td>
      <td><p align="center">0.75</p></td>
      <td><p align="center">4</p></td>
      <td><p align="center">1</p></td>
    </tr>
    <tr>
      <td><p align="center">5.0</p></td>
      <td><p align="center">0.75</p></td>
      <td><p align="center">4</p></td>
      <td><p align="center">1</p></td>
    </tr>
  </table>
</div>
<p align="center">
  <i>
    Table 5: Summary of algorithm performance based on gap penalty.
  </i>
</p>


## 🔗References

JONES, D. Threader: protein sequence threading by double dynamic programming. *Computational Methods in Molecular Biology.* Elsevier, 1996. v. 32, cap. 13, p. 312–338


## ✉️Contact

For questions or issues, please open an issue on GitHub or contact [nadiajuckova@gmail.com](mailto:nadiajuckova@gmail.com).
