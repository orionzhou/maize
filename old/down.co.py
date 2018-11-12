#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import math 
import numpy as np
import os.path as op
import argparse

from bs4 import BeautifulSoup
from urllib import urlretrieve
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.common.exceptions import WebDriverException
import time
import web 

def parse_co_homepage(html, fo):
    soup = BeautifulSoup(html, "lxml")
    prods = soup.findAll("div", itemprop = "itemListElement")
    print len(prods)

    fho = open(fo, "w")
    print >>fho, "pid\tcol\tname\twprice\tnprice\turl"
    ptn1 = re.compile("^(\d+)\s+([\w\/]+)$")
    for prod in prods:
        lst = ptn1.findall(prod['data-itemid'])
        pid, col = lst[0]
        
        xn = prod.find("div", class_ = "product-name").find("a", class_ = "name-link")
        name, url = xn.get('title'), xn['href']

        wprice, nprice = "", ""
        xp = prod.find("div", class_ = "product-price")
        xp1 = xp.find("span", class_ = "price-sales")
        if xp1 == None:
            xp1 = xp.find("span", class_ = "standardprice")
            xp2 = xp.find("span", class_ = "salesprice")
            nprice = xp2.text
        wprice = xp1.text
        wprice = re.sub(r'[^0-9.]', '', wprice)
        nprice = re.sub(r'[^0-9.]', '', nprice)
        print >>fho, "\t".join([pid, col, name, wprice, nprice, url])
    fho.close()
def parse_co_prodpage(mid, html, do):
    soup = BeautifulSoup(html, "lxml")
    prod = soup.find("div", class_ = "pdp_main_container")
    if prod== None:
        return 0
    ptn1 = re.compile("url\(\"(http\S+)\"\)")
    ptnp = re.compile("$([\d,\.]+)$")
    ptne = re.compile("changeColor\(\'\w+\',\'(\d+)\',.*colorChangeOmni\(\'([\w \/]+)\'")
    
    note, detl = "", ""
    xn = prod.find("p", itemprop = "description")
    if not xn == None:
        note = xn.text.replace("\n", "")
        note = note.encode('ascii', 'ignore').decode('ascii')
    xd = prod.find("div", class_ = "pdp_description_content jspScrollable pdp_description_tabs_2 hidden")
    if not xd == None:
        detl = xd.find("ul").text.replace("\n", ";")
        detl = detl.encode('ascii', 'ignore').decode('ascii')
  
    colNum, col = '', ''
    x1 = prod.find("ul", class_ = "color_group_list color_swatch")
    x2 = x1.find_all("li")
    for x in x2:
        if x.has_attr('class'):
            jsp = x.find('img')['onclick']
            lst = ptne.findall(jsp)
            colNum, col = lst[0]
    
    wprice, nprice = "", ""
    x1 = prod.find("div", id = "productPrice")
    if x1 == None:
        print mid, "no price info"
        return 0
    x2 = x1.find("span", class_ = "price") 
    if x2 == None:
        xa = x1.find("span", class_ = "was_price")
        if xa == None:
            print mid, "no price info"
        xb = x1.find("span", class_ = "now_price")
        wprice = xa.text.split("$")[1]
        nprice = xb.text.split("$")[1]
    else:
        wprice = x2.text.split("$")[1]
    wprice.replace(",", "")
    nprice.replace(",", "")

    imgUrls = []
    x1 = prod.find_all("div", class_ = "s7thumb")
    if len(x1) > 0:
        for x in x1:
            img_url = ptn1.findall(x['style'])[0]
            img_url = img_url.split("?")[0]
            imgUrls.append(img_url)
    else:
        x1 = prod.find_all("div", id = "productDetailFullSizeImage_flyout")
        x2 = x1[0].find_all("img")
        img_url = x2[0]['src'].split("?")[0]
        imgUrls.append(img_url)
    for i in range(0, len(imgUrls)):
        imgName = "%s-%s-%d.jpg" % (mid, colNum, i+1)
        img_url = "%s?scl=2.0566667" % imgUrls[i]
        urlretrieve(img_url, op.join(do, imgName))
    return wprice, nprice, col, colNum, note, detl
def get_co_info(fi, fo, do):
    os.system("rm %s/*" % do)
    wd = webdriver.Firefox()
    ptne = re.compile("changeColor\(\'\w+\',\'(\d+)\',.*colorChangeOmni\(\'([\w \/]+)\'")
    
    ary = np.genfromtxt(fi, names = True, dtype = None, delimiter = "\t")
    fho = open(fo, "w")
    print >>fho, "mid\twprice\tnprice\tcol\tcolNum\tnote\tdetl"
    i, cids = 1, set()
    for row in ary:
        pid, mid, wid, name, url = row
#        if not mid == '30S6GRUM2L': continue
        if url[0] == "/":
            url = url_base + url
        if not wid == "US_" + mid:
            print "inconsistent:", pid, mid, wid
        
        wd.get(url)
        html = wd.page_source.encode("utf-8")
        soup = BeautifulSoup(html, "lxml")
        
        colNum, col = '', ''
        x1 = soup.find("ul", class_ = "color_group_list color_swatch")
        x2 = x1.find_all("li")
        jsps = []
        for x in x2:
            jsp = x.find('img')['onclick']
            lst = ptne.findall(jsp)
            colNum, col = lst[0]
            cid = "%s-%s" % (mid, colNum)
            if cid in cids:
                continue
            else:
                cids.add(cid)
                jsps.append(jsp)
        for jsp in jsps:
            wd.execute_script(jsp)        
            time.sleep(2)
            html = wd.page_source.encode("utf-8")
            wprice, nprice, col, colNum, note, detl = parse_mk_prodpage(mid, html, do)
            print i, mid, wprice, nprice, col, colNum
            print >>fho, "\t".join([mid, wprice, nprice, col, colNum, note, detl])
            i += 1
    fho.close()
    wd.quit()
        
dirw = '/home/youngn/zhoup/Data/web/co'
if not os.path.exists(dirw): os.makedirs(dirw)
url = "https://www.coach.com/shop/women-handbags?viewAll=true"
url_base = "https://www.coach.com"

if __name__ == '__main__':
    os.chdir(dirw)
    #html = web.get_html_hard(url)
    #parse_co_homepage(html, "01.tbl")
    get_co_info("01.tbl", "04.tbl", "05.imgs")
