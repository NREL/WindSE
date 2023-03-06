#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name ="windse",
      description="Wind Systems Engineering",
      version="2021.08.01",
      author="Ryan King, Jeffery Allen, Ethan Young, John Jasa",
      author_email="ryan.king@nrel.gov, jeff.allen@nrel.gov",
      url="https://github.com/NREL/WindSE",
      include_package_data=True,
      packages=['windse','windse_driver'],
      entry_points={'console_scripts' : ['windse = windse_driver.driver:main'],},

)