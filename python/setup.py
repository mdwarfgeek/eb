#!/usr/bin/env python

import os
import distutils.core
import numpy.distutils.misc_util

inc = numpy.distutils.misc_util.get_numpy_include_dirs()
inc.append("../src")

mod = distutils.core.Extension("eb",
                               include_dirs=inc,
                               library_dirs=["../src"],
                               libraries=["eb", "m"],
                               sources=["wrap.c"])

distutils.core.setup(name="eb",
                     version="0.01",
                     description="EB light curve generator",
                     ext_modules=[mod])

