#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math 
import os.path as op

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch, cm
from reportlab.lib import utils

from reportlab.lib.pagesizes import letter
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import PageBreak

data=  [['00', '01', '02', '03', '04'],
    ['10', '11', '12', '13', '14'],
    ['20', '21', '22', '23', '24'],
    ['30', '31', '32', '33', '34']]
t=Table(data,5*[0.4*inch], 4*[0.4*inch])
t.setStyle(TableStyle([
    ('BOX', (0,1), (1,2), 0.25, colors.black),
    ('VALIGN',(2,1),(3,2),'TOP'),
]))


doc = SimpleDocTemplate('ex.pdf', pagesize = letter, 
        rightMargin = 72, leftMargin = 72,
        topMargin = 72, bottomMargin = 18)
doc.build([t])
