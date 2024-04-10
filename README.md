## Traffic Assignment (TA) Tutorial
The only MSA + FW + PBA code implementation you can find in Python. The implementation is based on the book [Transportation Network Analysis](https://sboyles.github.io/blubook.html) by Stephen D. Boyles, Nicholas E. Lownes, and Avinash Unnikrishnan. 

Please note that the code is for educational purposes only. The code is not optimized for large-scale networks.

If you find the code helpful, please consider giving it a star.


### 0. Prerequisite
- Install [Git](https://git-scm.com/downloads)
- Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

### 1. Clone the repository

```bash
git clone https://github.com/sguo28/Traffic_Assignment_Tutorial.git
cd Traffic_Assignment_Tutorial
```
### 2. Pacakges installation

```bash
conda create --name ta python=3.9
conda activate ta
pip install -r requirements.txt
```

```bash
# (requirements.txt)
networkx==3.1
pandas==2.1.4
numpy==1.26.3
scipy==1.11.4
matplotlib==3.8.0
```

#### Install Jupiter Notebook

```bash
conda install -c conda-forge jupyterlab
```

#### Run Jupyter Notebook

```bash
jupyter lab
```

### 3. Run the Code
```bash
TA_MSA.ipynb
TA_FW.ipynb
```

### Supported Functions
- [x] Frank-Wolfe Algorithm
- [x] Method of Successive Averages (MSA)
- [x] Path-based Algorithm (PBA)
- [x] User Equilibrium (UE)
- [ ] System Optimal (SO)


