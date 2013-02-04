GALATEAS - mosesdecoder
============

Moses, the machine translation system

This repository contains the additional features added to Moses during the Galateas project:
- discovery manager: the XML-RPC server can provide information on his parameters and state, like "spoken" languages, moses.ini parameters, etc. Useful for building large scale Moses clusters.
- upload/download: a user can upload files to a given directory on the Moses server (think aditional training data for example) or download files from the server (feature used by the reranking framework to access parameter files linked to a given Moses instance)
- nbest retrieval: the XML-RPC server is now able to retrieve all the nbest for a given translation (use parameter "nBestSize")
- features: the XML-RPC adds also the translation features (individual scores for all the computed features) to the translation result map.

The GALATEAS project offers digital content providers an innovative approach to understanding users' behaviour by analysing language-based information from transaction logs and facilitates the development of improved navigation and search technologies for multilingual content access.

More details on: http://www.galateas.eu/
