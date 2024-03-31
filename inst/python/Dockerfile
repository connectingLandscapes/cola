FROM python

EXPOSE 8501
WORKDIR /rb_mig
#install pip requirements

ADD requirements.txt /rb_mig/requirements.txt
RUN pip install -r /requirements.txt


ADD . /rb_mig

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "app.py"]