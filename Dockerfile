FROM python:3.7.4-buster

ENV PYTHONUNBUFFERED 1

WORKDIR /home

# Setup python
RUN pip install --upgrade pip
COPY requirements.txt .
RUN pip install -r requirements.txt

# Jupyter notebook
RUN jupyter contrib nbextension install
RUN jupyter nbextension enable varInspector/main
RUN pip install autopep8
RUN jupyter nbextension enable code_prettify/autopep8
RUN jupyter nbextension enable hinterland/hinterland
RUN jupyter nbextension enable collapsible_headings/main
RUN jupyter nbextension enable highlight_selected_word/main
RUN jupyter nbextension enable splitcell/splitcell
RUN jupyter nbextension enable toc2/main
RUN jupyter nbextension enable select_keymap/main
RUN jt -t grade3 -cursc o -cursw 5 -T -cellw 95%
COPY ./support-files/notebook.json /root/.jupyter/nbconfig/notebook.json

COPY . .

CMD echo "\n\033[1m\033[44m\033[33mSuccessfully built!\033[0m\n⮑  Now visit \033[4mhttp:localhost:8888\033[0m for simulations ⚛\n" && jupyter notebook --ip=0.0.0.0 --no-browser --allow-root --port 8888 --NotebookApp.token='' --NotebookApp.password=''
