FROM python:3.8
# re: mkdir, https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=863199#23
COPY . .

# WORKDIR /hailproj
RUN mkdir -p /usr/share/man/man1 && \
    apt-get update && apt-get install -y \
    openjdk-11-jre-headless \
    && rm -rf /var/lib/apt/lists/* && \
    pip3 install -r requirements.txt

RUN pylint .


ENTRYPOINT [ "python", "./main.py" ]
