#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name ="windse",
      description="Wind Systems Engineering",
      version="2018.1.0_0.4",
      author="Ryan King, Jeffery Allen",
      author_email="ryan.king@nrel.gov",
      url="https://github.com/WISDEM/WindSE",
      py_modules=['windse', 'windse_driver'],
      entry_points={'console_scripts' : ['windse_driver = windse.windse_driver:main'],},

)