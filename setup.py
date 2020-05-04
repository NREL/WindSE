#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name ="windse",
      description="Wind Systems Engineering",
      version="2018.1.0_0.8",
      author="Ryan King, Jeffery Allen",
      author_email="ryan.king@nrel.gov",
      url="https://github.com/NREL/WindSE",
      include_package_data=True,
      packages=['windse','windse_driver'],
      entry_points={'console_scripts' : ['windse = windse_driver.driver:main'],},

)