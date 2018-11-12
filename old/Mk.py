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

def make_ad(fi, di, fo, tax, exc):
    fhi = open(fi, "r")
    doc = SimpleDocTemplate(fo, pagesize = letter, 
            rightMargin = 72, leftMargin = 72,
            topMargin = 72, bottomMargin = 18)
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name = 'mytitle', fontName = 'Helvetica-Bold', fontSize = 18, leading = 25, spaceafter = 25))
    styles.add(ParagraphStyle(name = 'mytext', leading = 16, alignment = TA_LEFT))
    styles.add(ParagraphStyle(name = 'mydetail', leading = 14, alignment = TA_LEFT))
    
    story = []
    fhi.readline()
    for line in fhi:
        lst = line.strip("\n").split("\t")
        if len(lst) < 10: break
        mid, name, wprice, nprice, col, colCode, imgNum, url, note, detl = lst
        imgNum = int(imgNum)
        wprice = wprice.replace(",", "")
        nprice = nprice.replace(",", "")
        if nprice == '': continue
        wprice, nprice = float(wprice), float(nprice)
        if nprice < 150: continue
        uprice = (wprice * (1+tax/100) + 60 + 20) * exc
        rprice = (nprice * (1+tax/100) + 60 + 20) * exc
        uprice = str(int(uprice))
        rprice = str(int(rprice))

        title = "%s" % name
        story.append(Paragraph(title, styles["mytitle"]))

        keys = ['MFID', 'Color', 'MSRP (RMB)', 'Sale (RMB)']
        vals  = [mid, col, uprice, rprice]
        ptext = "<br />".join(["<font size=10>%s:</font><font size=12> %s</font>" % (keys[i], vals[i]) for i in range(0, len(keys))])
        ptext = "<br />" + ptext
        ltxt = Paragraph(ptext, styles['mytext'])
       
        lines = detl.split(";")
        detl = "<br />".join(["<font size = 10>%s</font>" % l for l in lines])
        rtxt = Paragraph(detl, styles['mydetail'])

        t1 = Table([[ltxt, rtxt]])
        t1.setStyle(TableStyle([
            ('ALIGN', (0,0), (-1,-1), 'LEFT'),
            ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ]))
        story.append(t1)
        story.append(Spacer(1,12))

        imglayout = [['',''],['','']]
        for i in range(0, imgNum):
            imgpath = op.join(di, "%s-%s-%d.jpg" % (mid, colCode, i+1))
            iw, ih = utils.ImageReader(imgpath).getSize()
            aspect = ih / float(iw)
            wd = 3 
            ht = wd * aspect
            img = Image(imgpath)#, wd*inch, ht*inch)
            if i == 0:
                imglayout[0][0] = img
            elif i == 1:
                imglayout[0][1] = img
            elif i == 2:
                imglayout[1][0] = img
            else:
                imglayout[1][1] = img
        t = Table(imglayout)
        story.append(t)
        story.append(PageBreak())
    doc.build(story)
    return True
def easy_draw():
    dirw = '/home/youngn/zhoup/Data/web'
    fignames = ["%d.jpg" % num for num in range(1,5)]
    figs = [op.join(dirw, figname) for figname in fignames]
    
    c = canvas.Canvas(op.join(dirw, 'ex2.pdf'))
    for i in range(0,4):
        img = utils.ImageReader(figs[i])
        iw, ih = img.getSize()
        aspect = ih / float(iw)
        wd = 10 
        ht = wd * aspect
        gap = 0.5
        x = (wd+gap)*i*cm
#        c.drawImage(figs[i], x, wd*cm, width = wd*cm, height = wd*aspect*cm)
        if i == 0:
            c.drawImage(figs[i], gap*2*cm, (gap*2+wd)*cm, width = wd*cm, height = ht*cm)
        elif i == 1:
            c.drawImage(figs[i], (gap*2+wd)*cm, (gap*2+wd)*cm, width = wd*cm, height = ht*cm)
        elif i == 2:
            c.drawImage(figs[i], gap, gap, width = wd*cm, height = ht*cm)
        elif i == 3:
            c.drawImage(figs[i], (gap*2+wd)*cm, gap, width = wd*cm, height = ht*cm)
        else:
            print "error: %d figs" % len(figs)
    c.showPage()
    c.save()

if __name__ == '__main__':
    dirw = '/home/youngn/zhoup/Data/web/mk'
    os.chdir(dirw)
    tax = 7.625
    exc = 6.6
    make_ad("04.tbl", "05.imgs", "10.pdf", tax, exc)
