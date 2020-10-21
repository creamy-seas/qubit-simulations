# Setup #

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
