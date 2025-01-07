#!/usr/bin/env python3
import os
import subprocess
from time import sleep
cwd = os.getcwd() #current working directory
lwd = os.getcwd() + "/lib" #library (Annotation_libraries) working directory
twd = os.getcwd() + "/bin" #tools (Assembly_Tools) working directory
ori_tmm = '/xxxxx'
#skipped:python3 installation, Oases installation, Shannon installation

def main():
    subprocess.run(f"sed -i 's!{ori_tmm}!{cwd}!g' {twd}/busco/config/config.ini", shell=True)
    print("\n\nBusco's config.ini File Correctly Configured!")
    sleep(3)
    subprocess.run(f"sudo apt install curl", shell=True)
    subprocess.run(f"curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py", shell=True)
    subprocess.run(f"python get-pip.py", shell=True)
    subprocess.run(f"sudo apt install python3-pip", shell=True)
    #subprocess.run(f"sudo apt install python-pip", shell=True)    
    subprocess.run(f"sudo pip3 install numpy", shell=True)
    subprocess.run(f"sudo pip install numpy", shell=True)
    subprocess.run(f"sudo apt install gcc", shell=True)
    subprocess.run(f"sudo apt install make", shell=True)
    subprocess.run(f"sudo apt install cmake", shell=True)
    subprocess.run(f"sudo apt install libbz2-dev", shell=True)
    subprocess.run(f"sudo apt install libncurses5-dev", shell=True)
    subprocess.run(f"sudo apt install zlib1g-dev", shell=True)
    subprocess.run(f"sudo apt install libncursesw5-dev", shell=True)
    subprocess.run(f"sudo apt install liblzma-dev", shell=True)
    subprocess.run(f"sudo apt install libcurl4-openssl-dev", shell=True)
    subprocess.run(f"sudo apt install libssl-dev", shell=True)
    subprocess.run(f"sudo apt install libsparsehash-dev", shell=True)
    subprocess.run(f"sudo pip3 install shmlast", shell=True)
    subprocess.run(f"sudo apt install bowtie2", shell=True)
    subprocess.run(f"sudo apt install samtools", shell=True)
    subprocess.run(f"pip3 install pandas", shell=True)
    subprocess.run(f"pip3 install biopython", shell=True)
    subprocess.run(f"pip3 install tqdm", shell=True)
    #subprocess.run(f"sudo apt install sra-toolkit", shell=True)
    subprocess.run(f"sudo apt-get update -y", shell=True)
    subprocess.run(f"sudo apt install default-jre", shell=True)
    subprocess.run(f"sudo apt install git", shell=True)
    subprocess.run(f"sudo apt install dos2unix", shell=True)
    subprocess.run(f"sudo apt install cd-hit", shell=True)
    subprocess.run(f"sudo apt install pkgconf", shell=True)
    subprocess.run(f"sudo pip3 install sh", shell=True)
    subprocess.run(f"sudo pip install cvxopt", shell=True)
    subprocess.run(f"sudo pip3 install cvxopt", shell=True)
    subprocess.run(f"sudo apt install metis", shell=True)
    subprocess.run(f"sudo apt-get -y install r-base", shell=True)
    subprocess.run(f"pip install python-igraph", shell=True)
    subprocess.run(f"sudo apt-get install csvtool", shell=True)
    os.chdir(twd+"/Rcorrector-1.0.4")
    subprocess.run(f"make", shell=True)
    os.chdir(twd+"/velvet_1.2.10")
    subprocess.run(f"make 'MAXKMERLENGTH=149' 'CATEGORIES=5'", shell=True)
    os.chdir(twd+"/oases-0.2.09/oases")
    velvet = f"{twd}/velvet_1.2.10"
    subprocess.run(f"make 'VELVET_DIR={velvet}' 'MAXKMERLENGTH=149' 'CATEGORIES=5'", shell=True)
    os.chdir(twd+"/abyss-2.2.4")
    boost = f"{twd}/boost_1_72_0"
    #mpi = f"{twd}/openmpi-4.0.4"
    subprocess.run(f"./configure --with-boost={boost}", shell=True)
    #subprocess.run(f"./configure --with-mpi={mpi}", shell=True)
    subprocess.run(f"make", shell=True)
    subprocess.run(f"sudo make install", shell=True)
    os.chdir(twd+"/jellyfish-2.3.0")
    subprocess.run(f"./configure", shell=True)
    subprocess.run(f"make", shell=True)
    subprocess.run(f"sudo make install", shell=True)
    os.chdir(twd+"/trinityrnaseq-v2.11.0.FULL/trinityrnaseq-v2.11.0")
    subprocess.run(f"make", shell=True)
    subprocess.run(f"make plugins", shell=True)
    subprocess.run(f"sudo make install", shell=True)
    os.chdir(twd+"/ncbi-blast-2.10.0+-x64-linux/ncbi-blast-2.10.0+")
    subprocess.run(f"export PATH=$PATH:$HOME/{twd}/ncbi-blast-2.10.0+-x64-linux/ncbi-blast-2.10.0+", shell=True) #append export location
    os.chdir(twd+"/RSEM-1.3.1")
    subprocess.run(f"make", shell=True)
    subprocess.run(f"sudo make install", shell=True)
    os.chdir(twd+"/busco")
    subprocess.run(f"python3 setup.py install --user", shell=True)
    os.chdir(twd+"/hmmer/hmmer-3.2.1")
    subprocess.run(f"./configure", shell=True)
    subprocess.run(f"make", shell=True)
    subprocess.run(f"sudo make install", shell=True)
    os.chdir(twd+"/TransLiG_1.3")
    boost = f"{twd}/boost_1_72_0"
    subprocess.run(f"./configure --with-boost={boost}", shell=True)
    subprocess.run(f"make", shell=True)
    os.chdir(twd+"/BinPacker_1.0")
    boost = f"{twd}/boost_1_72_0"
    subprocess.run(f"./configure --with-boost={boost}", shell=True)
    subprocess.run(f"make", shell=True)
    os.chdir(cwd)
    subprocess.run(f"sudo apt-get -y install libboost-iostreams-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install libgsl-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install libboost-graph-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install libboost-all-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install libsuitesparse-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install liblpsolve55-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install libsqlite3-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install libmysql++-dev", shell=True)
    subprocess.run(f"sudo apt-get -y install libbamtools-dev", shell=True)
    os.chdir(twd+"/Augustus-master")
    subprocess.run(f"make", shell=True)
    subprocess.run(f"sudo make install", shell=True)
    os.chdir(cwd)
    subprocess.run(f"sudo apt install python", shell=True)
    subprocess.run(f"sudo apt-get install -y pigz", shell=True)
    subprocess.run(f"sudo apt-get install python3-tk", shell=True)
    #subprocess.run(f"sudo pip install python-igraph", shell=True)#pip not working, use conda?
    subprocess.run(f"sudo pip3 install python-igraph", shell=True)
    os.chdir(twd)
    setup_apt = f"{twd}/setup-apt.sh"
    subprocess.run(f"sudo {setup_apt}", shell=True)
    subprocess.run(f"source /etc/profile.d/sra-tools.sh", shell=True)
    subprocess.run(f"vdb-config --interactive", shell=True)
    os.chdir(cwd)
    subprocess.run(f"sudo ldconfig", shell=True)
    print("\n\nInstall Completed\n\n")
    sleep(3)
    
if __name__ == '__main__':
    main()


######################################################################################################################

######################################################################################################################
                                                ###   END April 2 2021
######################################################################################################################

######################################################################################################################
