FROM python:3.7-slim as base

# Setup environment
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONFAULTHANDLER 1

FROM base AS deps

# Install pipenv and gcc
RUN pip install pipenv
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc

# Install python dependencies in /.venv
COPY Pipfile .
COPY Pipfile.lock .
RUN PIPENV_VENV_IN_PROJECT=1 pipenv install --deploy

FROM base AS runtime

# Copy virtual env from python-deps stage
COPY --from=deps /.venv /.venv
ENV HDF5_PLUGIN_PATH="example/deps/"
ENV PATH="/.venv/bin:$PATH"

# Install application into container
RUN mkdir -p /app
WORKDIR /app
COPY . .
