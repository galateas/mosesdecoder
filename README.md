GALATEAS - mosesdecoder
============

Moses, the machine translation system

This repository contains features added to Moses in the context of the Galateas project:
- discovery manager: the XML-RPC server can provide information on its parameters and state, like list of languages supported, moses.ini parameters, etc. Useful for building large scale Moses clusters.
- upload/download: a user can upload files to a given directory on the Moses server (such as additional training data for example) or download files from the server.
- nbest retrieval: the XML-RPC server is now able to retrieve all the nbest translation hypotheses for a given translation (use parameter "nBestSize").
- translation parameters: the XML-RPC adds also the translation parameters (individual scores for all the computed features) to the translation result map.

The GALATEAS project
============
The GALATEAS project offers digital content providers an innovative approach to understanding users' behaviour by analysing language-based information from transaction logs and facilitates the development of improved navigation and search technologies for multilingual content access.

More details on: http://www.galateas.eu/
