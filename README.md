# Setup Manual #
> Will need to install python and all dependecies manually

1. Create python virtual environment
2. Install python modules
```shell
pip install -r requirements.txt
```

3. Install juputer modules
```shell
jupyter contrib nbextension install --user
jupyter nbextension enable hinterland/main
jupyter nbextension enable varInspector/main
jupyter nbextension enable code_prettify/autopep8
jupyter nbextension enable collapsible_headings/main
jupyter nbextension enable toc2/main
jupyter nbextension enable highlight_selected_word/main
jupyter nbextension enable splitcell/splitcell
jupyter nbextension enable select_keymap/main
jupyter nbextension enable toggle_all_line_numbers/main
jupyter nbextension enable autosavetime/main

jt -t grade3 -cellw 70%  -altm -cursc g -cursw 5 -T
```

4. Run project
```shell
jupyter notebook
```


# Setup Automated #
> - Will work on all machines
> - Will need to install `docker`

## Installation ##

> See [this video](https://drive.google.com/file/d/1U5oR8yTkWLR3nNUIOpEvF3HtEKRtJ8bG/view?usp=sharing) oh how to install.

1. Install [**Docker desktop**](https://docs.docker.com/desktop/)

2. Launch docker application

3. Clone this directory - in the command line write
```shell
git clone https://github.com/creamy-seas/qubit-simulations
```

4. Go to this directory and start the program. It will take 5 minutes to build.
```
cd qubit-simulations

make
```

5. Type in http://localhost:8888 in your browser

## Running ##

After running the above installation, please see [this video](https://drive.google.com/file/d/1ia9tHd4D7tmGBfza3BAfx8aNT4bh4AsL/view?usp=sharing) on how to start and stop the program.

# Examples #
<details>
<summary>Shapiro simulations</summary>

> File: [2021-03_shapiro-step-simulations.ipynb](./2021-03_shapiro-step-simulations.ipynb)

![shapiro-v1](./support-files/2021-03-11(Thu)_shapiro-simulation-v1.gif)

![shapiro-v2](./support-files/2021-03-11(Thu)_shapiro-simulation-v2.gif)


</details>
