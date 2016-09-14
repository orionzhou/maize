#!/usr/bin/env python
# -*- coding: utf-8 -*-
from BeautifulSoup import BeautifulSoup
import urllib2
import re

url = "http://www.loldytt.com/Anime/haizeiwang/"
url = "http://www.loldytt.com/Anime/HZWJCBHJ/"
if __name__ == '__main__':
    html_page = urllib2.urlopen(url)
    soup = BeautifulSoup(html_page)
    for link in soup.findAll('a'):
        print link.get('href')
