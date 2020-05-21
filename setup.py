#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name ="windse",
      description="Wind Systems Engineering",
      version="2020.5.0",
      author="Ryan King, Jeffery Allen, Ethan Young",
      author_email="ryan.king@nrel.gov, jeff.allen@nrel.gov",
      url="https://github.com/NREL/WindSE",
      include_package_data=True,
      packages=['windse','windse_driver'],
      entry_points={'console_scripts' : ['windse = windse_driver.driver:main'],},

)