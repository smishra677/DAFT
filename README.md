
# DAFT Tutorial

## DAFT: Discordant Attachment Frequency Test

DAFT (Discordant Attachment Frequency Test) is designed to distinguish
**introgression** from **incomplete lineage sorting (ILS)** using gene trees
and a species tree.

---

## DAFT Significance
**Purpose:** Distinguishing introgression from ILS.

### Step 1: Results without Loss

- **Script:** `Significance.py`

#### Input
- **SPECIES_TREE (`sp`)**  
  Species tree corresponding to the set of gene trees.

- **LIST_OF_GENE_TREES (`gt`)**  
  List of gene trees to test.

- **Lineage1/Lineage2 (`lineage`)**  
  Two lineages to test lineage loss, separated by `/`.  
  Can also be set to `NONE`.

- **sibling**  
  Flag to perform sibling test:  
  `1` = perform sibling test  
  `0` = do not perform sibling test

#### Output
- **`DAFT_Significance.txt`**  
  (See **DAFT_Significance_Explanation** below.)

#### Command
```bash
python Significance.py --sp <SPECIES_TREE> --gt <LIST_OF_GENE_TREES> --lineage NONE --sibling <1/0>
````

---

### Step 2: Results with Loss

* **Script:** `Significance.py` (same as Step 1)

* **Lineage1/Lineage2**
  Provide two lineages (not `NONE`), for example: `1;/2;`

#### Output

* **`DAFT_Significance.txt`**
  (See **DAFT_Significance_Explanation** below.)

#### Command

```bash
python Significance.py --sp <SPECIES_TREE> --gt <LIST_OF_GENE_TREES> --lineage <Lineage1/Lineage2> --sibling <1/0>
```

---

## DAFT Direction

**Purpose:** Detecting receiver and donor lineages.

### Step 1

* **Script:** `Direction.py`

#### Input

* **SPECIES_TREE (`sp`)**
  Species tree corresponding to the set of gene trees.

* **LIST_OF_GENE_TREES (`gt`)**
  List of gene trees to test.

* **SIGNIFICANT_PAIRS (`sig`)**
  List of lineage pairs identified as significant.

#### Output

* **`DAFT_Direction.txt`**
  (See **DAFT_Direction_Explanation** below.)

#### Command

```bash
python Direction.py --sp <SPECIES_TREE> --gt <LIST_OF_GENE_TREES> --sig <LIST_OF_SIGNIFICANT_PAIRS>
```

---

## DAFT_Significance_Explanation

### Snapshot of `DAFT_Significance.txt`

```
    Species Tree = ((((1,2),3),((4,5),6)),7);
    ========================================
    Focal_lineage = (2,1);
    NNI_sp  Test_lineage  total_count                   comparison_uncle  uncle_count  Z-value-uncle                    comparison_sibling  sibling_count  Z-value-sibling
    ==================================================================================================================================================================================================
    0       3;          1754                                          -            -                                                    -              -              
    3       4;          11                          6;                0            -3.32                            5;                  0              -3.32          
    --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

```

### Field Descriptions

* **Species Tree**
  Species tree provided as input.

* **Focal_lineage**
  Branch where counts are evaluated.

* **Test_lineage**
  Lineage counted with respect to the focal lineage.

* **NNI_sp**
  Number of NNIs separating the test lineage from the focal lineage
  in the species tree.

* **total_count**
  Number of times `(Focal_lineage, Test_lineage)` occur together
  in the gene trees.

* **comparison_uncle**
  Avuncular lineage of the test lineage in the species tree.

* **uncle_count**
  Number of times `(Focal_lineage, comparison_uncle)` occur together.

* **Z-value-uncle**
  Z-score comparing `uncle_count` and `total_count`.

* **comparison_sibling**
  Sibling lineage of the test lineage (if sibling test is enabled).

* **sibling_count**
  Number of times `(Focal_lineage, comparison_sibling)` occur together.

* **Z-value-sibling**
  Z-score comparing `sibling_count` and `total_count`.

More negative Z-values indicate stronger statistical significance.

---

## DAFT_Direction_Explanation


```
Significant_Pairs
=========================================
Between 1; AND 4;           4  NNI APART                 
Between 4; AND (2,1);       3  NNI APART 
Between 1; AND (5,4);       3  NNI APART  
Between 1; AND 5;           4  NNI APART          
Between 4; AND 2;           4  NNI APART 
*****************************************


DATA TABLE1
=====================================================================================
Significant_Pairs  Total_gene_trees  Lineage1  Count1  Lineage2  Count2  Major_Moved 
=====================================================================================
1; AND 4;          129               1;        46      4;        83      4;   		                  
4; AND (2,1);      11                4;        11      (2,1);    0       4;          
1; AND (5,4);      26                1;        26      (5,4);    0       1;            
1; AND 5;          31                1;        31      5;        0       1;              
4; AND 2;          11                4;        11      2;        0       4;           
*************************************************************************************



DATA TABLE2(BIDIRECTIONAL)
=========================================================================================================================
Significant_Pairs  Minor_Moved   Minor_sibling  Total_gene_trees    Minor_sibling_count   Minor_moved_count     Z_score
=========================================================================================================================
1; AND 4;          1;		     2;	            11        			0                     46     		       -6.78                   
*************************************************************************************************************************


Inferred_Relations:
==========================================================================
Between 1; AND 4;        Receiver: 4; AND Donor: 1;         (BIDIRECTIONAL)          
Between 4; AND (2,1);    Receiver: 4; AND Donor: (2,1);
Between 1; AND (5,4);    Receiver: 1; AND Donor: (5,4);
Between 1; AND 5;        Receiver: 1; AND Donor: 5;      
Between 4; AND 2;        Receiver: 4; AND Donor: 2;
***************************************************************************

==========================================================================================================
Network: (7,((6,(#H2,((#H5,5),((((4)#H1)#H3)#H6,#H4)))),(3,(#H1,((#H6,2),((((1)#H2)#H4)#H5,#H3))))));
********************************************************************************************************** 
```


### Significant Pairs

Significant Pairs that were passed as input:
"Between 1; AND 4;           4  NNI APART" means lineage 1; and lineage 4; were pass as input and are 4 NNIs apart in species tree.


```
Between 1; AND 4;           4 NNI APART
Between 4; AND (2,1);       3 NNI APART
Between 1; AND (5,4);       3 NNI APART
Between 1; AND 5;           4 NNI APART
Between 4; AND 2;           4 NNI APART
```

---

### Data Table 1: Unidirectional Introgression Test

**Columns**

* **Significant_Pairs**: Input lineage pair.
* **Total_gene_trees**: Number of gene trees containing the pair.
* **Lineage1 / Lineage2**: Lineages in the pair.
* **Count1 / Count2**: Number of times each lineage dominates NNI transformations.
* **Major_Moved**: Lineage with the higher count (receiver lineage).

---

### Data Table 2: Bidirectional Test

* **Minor_Moved**
  Lineage with the lower count from Data Table 1.

* **Minor_sibling**
  Sibling of `Minor_Moved` in the species tree.

* **Total_gene_trees**
  Number of trees containing `(Major_Moved, Minor_sibling)`.

* **Minor_sibling_count**
  Number of times `Minor_sibling` dominates.

* **Minor_moved_count**
  Count of `Minor_Moved`.

* **Z_score**
  Z-score comparing `Minor_sibling_count` and `Minor_moved_count`.

More negative Z-scores indicate stronger statistical significance.

---

## Inferred Relations

Each significant pair is assigned a **Receiver** and **Donor** lineage.
Pairs are marked as **BIDIRECTIONAL** if supported by bidirectional tests.

---

## Network

The inferred introgression network is represented as a Newick string
and can be visualized using tools such as **ICYTREE** or **PHYLOPLOTS**.

```
(7,((6,(#H2,((#H5,5),((((4)#H1)#H3)#H6,#H4)))),(3,(#H1,((#H6,2),((((1)#H2)#H4)#H5,#H3))))));
```

```
