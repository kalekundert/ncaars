FROM python:3.10.2-alpine3.15

# NextFlow seems to require bash
RUN apk add --no-cache bash

COPY rosetta /rosetta
ENV PATH="/rosetta/source/scripts/python/public:${PATH}"
