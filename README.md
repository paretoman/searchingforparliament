# Searching for Parliament

1. Open a terminal or cmd to where you want to work.

2. Download this project.
    ```
	git clone https://github.com/paretoman/searchingforparliament.git
	```

3. Install Gurobi (I only got it working with a [Miniconda](https://conda.io/miniconda.html) install of python2)
   * Sign up for an account for free [link](http://www.gurobi.com/registration/download-reg)
   *  Log in
   * [Download gurobi](http://www.gurobi.com/downloads/download-center)
    ```
	conda config --add channels http://conda.anaconda.org/gurobi
	conda install gurobi
    ```
   * [Request license](https://user.gurobi.com/download/licenses/free-academic)	
   * Go to your terminal and do the grbgetkey command that you got from that website.  When it asks where to save the key, save it in the default location. ( there is also a [quick start guide](https://www.gurobi.com/documentation/7.0/quickstart_windows/index.html) that has these steps.  And for [linux](https://www.gurobi.com/documentation/7.0/quickstart_linux/index.html) and [mac](https://www.gurobi.com/documentation/7.0/quickstart_mac/index.html))
   * You have to be physically on a campus for gurobi to work because it checks your ip.  Or you can use something like [sshuttle](https://superuser.com/questions/62303/how-can-i-tunnel-all-of-my-network-traffic-through-ssh) to forward all your traffic to the campus if you can log in.

4. Install openstv (my custom version on [github](https://github.com/paretoman/openstv) forked from [old](https://github.com/Conservatory/openstv) ).  My custom version makes it so we don't have to write and read files to use openstv.
    ```
	git clone https://github.com/paretoman/openstv.git
	cd openstv
    python setup.py install
    ```

5. Start Python's webserver from the command line

```
make
```

or

```
python SimpleServer.py
```

or on windows

```
c:\Users\You\Miniconda2\python.exe SimpleServer.py
```

6. Point your browser at http://localhost:8000

