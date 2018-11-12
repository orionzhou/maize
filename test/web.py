#!/usr/bin/env python
# -*- coding: utf-8 -*-
from urllib2 import urlopen
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.common.exceptions import WebDriverException
import time

def internet_on():
    try:
        resp = urllib2.urlopen("http://www.google.com/", timeout = 3)
        return True
    except urllib2.URLError as err: pass
    return False
def ajax_complete(driver):
    try:
        return 0 == driver.execute_script("return jQuery.active")
    except WebDriverException:
        pass
def get_html_easy(url):
    return urlopen(url).read()
def get_html_medium(url):
    wd = webdriver.Firefox()
    wd.get(url)
    html = wd.page_source.encode("utf-8")
    wd.quit()
    return html
def get_html_hard(url):
    #wd = webdriver.Firefox()
    wd = webdriver.PhantomJS()
    wd.set_window_size(800, 600)
    wd.get(url)
    
    lastHeight = wd.execute_script("return document.body.scrollHeight")
    rd = 1
    print rd, lastHeight
    while True:
        wd.execute_script("window.scrollTo(0, document.body.scrollHeight);")    
        time.sleep(3)
        WebDriverWait(wd, 60).until(ajax_complete, "Timeout waiting")
        newHeight = wd.execute_script("return document.body.scrollHeight")
        if newHeight == lastHeight: break
        lastHeight = newHeight
        rd += 1
        print rd, lastHeight
    html = wd.page_source.encode("utf-8")
    wd.quit()
    return html
       
if __name__ == '__main__':
    print "You lost your way"
