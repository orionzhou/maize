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
import Mk

def parse_mk_homepage(html, fo):
    soup = BeautifulSoup(html, "lxml")
    prods = soup.findAll("li", "product_panel_medium product_panel_data")
    print len(prods)

    fho = open(fo, "w")
    print >>fho, "pid\tmid\twid\tname\turl"
    ptn1 = re.compile("^[^$]*\$([\d,\.]+)$")
    ptn2 = re.compile("\$([\d,\.]+)[\n\r]*\-[\n\r]*\$([\d,\.]+)$")
    for prod in prods:
        pid, mid, wid = prod['data-productid'], prod['data-mfritemnum'], prod['id']
        name = prod.find("h6").text.encode('ascii', 'ignore').decode('ascii')
        wprice, nprice = "", ""
        xp = prod.find("span", "price")
        wprice.replace(",", "")
        nprice.replace(",", "")
        x1 = prod.find("div", "product_panel")
        urls = x1.find_all('a', oncontextmenu="javascript:mkobj.setCurrentURLInCookie(this);")
        url = urls[0].get('href').encode('ascii', 'ignore').decode('ascii')
        print >>fho, "\t".join([pid, mid, wid, name, url])
    fho.close()
def parse_mk_prodpage(mid, soup, do):
    prod = soup.find("div", class_ = "pdp_main_container")
    if prod == None: return 0
    ptn1 = re.compile("background-image: url\(\"?(http\S+)\"?\)")
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
  
    col, colCode = '', ''
    x1 = prod.find("ul", class_ = "color_group_list color_swatch")
    x2 = x1.find_all("li")
    for x in x2:
        if x.has_attr('class'):
            jsp = x.find('img')['onclick']
            lst = ptne.findall(jsp)
            colCode, col = lst[0]
    
    wprice, nprice = "", ""
    x1 = prod.find("div", id = "productPrice")
    if x1 == None: return False, []
    x2 = x1.find("span", class_ = "price") 
    if x2 == None:
        xa = x1.find("span", class_ = "was_price")
        if xa == None: return False, []
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
            lst = ptn1.findall(x.get('style'))
            if len(lst) == 0: return False, []
            img_url = lst[0].split("?")[0]
            imgUrls.append(img_url)
    else:
        x1 = prod.find("div", id = "productDetailFullSizeImage_flyout")
        x2 = x1.find_all("img")
        if not x2[0].has_attr('src'): return False, []
        img_url = x2[0]['src'].split("?")[0]
        imgUrls.append(img_url)
    for i in range(0, len(imgUrls)):
        imgName = "%s-%s-%d.jpg" % (mid, colCode, i+1)
        img_url = "%s?scl=2.0566667" % imgUrls[i]
        urlretrieve(img_url, op.join(do, imgName))
    return wprice, nprice, col, colCode, len(imgUrls), note, detl

def check_mkprod_html(html):
    soup = BeautifulSoup(html, "lxml")
    ptn1 = re.compile("background-image: url\(\"?(http\S+)\"?\)")
    ptne = re.compile("changeColor\(\'\w+\',\'(\d+)\',.*colorChangeOmni\(\'([\w \/]+)\'")
    
    prod = soup.find("div", class_ = "pdp_main_container")
    if prod == None: 
        main = soup.find("div", id = "main_container")
        if not main == None:
            if "product is no longer available" in main.text:
                return 2, soup
        return 1, soup 
    
    x1 = prod.find("div", id = "productPrice")
    if x1 == None:
        return 1, soup
    x2 = x1.find("span", class_ = "price") 
    if x2 == None:
        xa = x1.find("span", class_ = "was_price")
        if xa == None:
            return 1, soup

    x1 = prod.find_all("div", class_ = "s7thumb")
    if len(x1) > 0:
        for x in x1:
            lst = ptn1.findall(x.get('style'))
            if len(lst) == 0:
                return 1, soup
    else:
        x1 = prod.find("div", id = "productDetailFullSizeImage_flyout")
        if x1 == None:
            return 1, soup
        x2 = x1.find_all("img")
        if not x2[0].has_attr('src'):
            return 1, soup 
    return 0, soup
def get_soup_mkprod(wd, url, jsp, mid):
    flag = 1
    while flag == 1:
        wd.get(url)
        WebDriverWait(wd, 60).until(web.ajax_complete, "Timeout waiting")
        html = wd.page_source.encode("utf-8")
        flag, soup = check_mkprod_html(html)
        if flag == 1:
            print "\terror requesting", mid
            continue
        elif flag == 2:
            print "\tno longer available", mid
        if not jsp == None:
            wd.execute_script(jsp)
            WebDriverWait(wd, 60).until(web.ajax_complete, "Timeout waiting")
            html = wd.page_source.encode("utf-8")
        flag, soup = check_mkprod_html(html)
        if flag == 1:
            print "\terror requesting", mid
    return flag, soup
def pull_mk_info(fi, fo, do):
    os.system("rm %s/*" % do)
    wd = webdriver.Firefox()
    #wd = webdriver.PhantomJS()
    wd.set_window_size(800, 600)
    ptne = re.compile("changeColor\(\'\w+\',\'(\d+)\',.*colorChangeOmni\(\'([\w \/]+)\'")
    
    ary = np.genfromtxt(fi, names = True, dtype = None, delimiter = "\t")
    fho = open(fo, "w")
    print >>fho, "mid\tname\twprice\tnprice\tcol\tcolCode\timgNum\turl\tnote\tdetl"
    i, cids = 1, set()
    for row in ary:
        pid, mid, wid, name, url = row
#        if not mid == '30S4GTVT2L': continue
        if url[0] == "/":
            url = url_base + url
        if not wid == "US_" + mid:
            print "inconsistent:", pid, mid, wid
        flag, soup = get_soup_mkprod(wd, url, None, mid)
        if flag == 2: continue
        
        colNum, col = '', ''
        x1 = soup.find("ul", class_ = "color_group_list color_swatch")
        if x1 == None: continue
        x2 = x1.find_all("li")
        opts = []
        for x in x2:
            jsp = x.find('img')['onclick']
            lst = ptne.findall(jsp)
            colNum, col = lst[0]
            
            if x.has_attr("class"):
                wprice, nprice, col, colCode, imgNum, note, detl = parse_mk_prodpage(mid, soup, do)
                print i, mid, wprice, nprice, col, colNum, imgNum
                print >>fho, "\t".join([mid, name, wprice, nprice, col, colCode, str(imgNum), url, note, detl])
                i += 1
            else:
                opts.append([colNum, col, jsp])

            cid = "%s-%s" % (mid, colNum)
            if cid in cids:
                continue
            else:
                cids.add(cid)
        for opt in opts:
            colNum, col, jsp = opt
            flag, soup = get_soup_mkprod(wd, url, jsp, mid)
            if flag == 2: continue
            wprice, nprice, col, colCode, imgNum, note, detl = parse_mk_prodpage(mid, soup, do)
            print i, mid, wprice, nprice, col, colNum, imgNum
            print >>fho, "\t".join([mid, name, wprice, nprice, col, colCode, str(imgNum), url, note, detl])
            i += 1
    fho.close()
    wd.quit()
        
dirw = '/home/youngn/zhoup/Data/web/mk'
if not os.path.exists(dirw): os.makedirs(dirw)
url = "http://www.michaelkors.com/women/handbags/_/N-28f3"
url = "http://www.michaelkors.com/handbags/view-all-handbags/_/N-283i"
url_base = "http://www.michaelkors.com"

if __name__ == '__main__':
    os.chdir(dirw)
    #html = web.get_html_hard(url)
    #parse_mk_homepage(html, "01.tbl")
    #pull_mk_info("01.tbl", "04.tbl", "05.imgs")
    tax = 7.625
    exc = 6.6
    Mk.make_ad("04.tbl", "06.imgs.wm", "10.pdf", tax, exc)
    
