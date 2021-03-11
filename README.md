# Setup Manual #
> Will need to install python and all dependecies manually.
> See [**this video**](https://drive.google.com/file/d/1W-4pjqTsHDLQ3LV9SSvA2kXKPJqugf_G/view?usp=sharing) for example setup.

1. Clone this directory - in the command line write
```shell
git clone https://github.com/creamy-seas/qubit-simulations
```
2. Create python virtual environment
3. Install python modules
```shell
pip install -r requirements.txt
```

4. Install juputer modules
```shell
jupyter contrib nbextension install
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

5. Run project
```shell
jupyter notebook
```

# Setup that will work on all computers #

## Installation (only once) ##

> See [this video](https://drive.google.com/file/d/1ZFS4tugP9XYUNOpvOWhaisQyw6EHkO9B/view?usp=sharing) oh how to install.

1. Install [**Docker desktop**](https://docs.docker.com/desktop/)

2. Launch docker application

3. Clone this directory - in the command line write
```shell
git clone https://github.com/creamy-seas/qubit-simulations
```

4. Go to this directory and start the program. It will take 5 minutes to build.
```shell
cd qubit-simulations

make
```

## Running ##

1. Go to this directory and run (if you just ran the installation, then program is already running and this is not needed)

```shell
make
```

2. Type in http://localhost:8888 in your browser to open the project

Alternatively, see how the program is launched in [this video](https://drive.google.com/file/d/1yEXCtdDS1q6IzbYI_V0LQ01JDc0eh_pO/view?usp=sharing).

# Examples #
<details>
<summary>Shapiro simulations</summary>

> File: [2021-03_shapiro-step-simulations.ipynb](./2021-03_shapiro-step-simulations.ipynb)

![shapiro-v1](./support-files/2021-03-11(Thu)_shapiro-simulation-v1.gif)

![shapiro-v2](./support-files/2021-03-11(Thu)_shapiro-simulation-v2.gif)

![shapiro-step](./support-files/2021-03-11(Thu)_shapiro-simulation-v3.png)


</details>
