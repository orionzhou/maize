#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials

scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']
f_cred = '/home/springer/zhoux379/.config/google_account_token.json'
cred = ServiceAccountCredentials.from_json_keyfile_name(f_cred, scope)
gc = gspread.authorize(cred)

#book = gc.open_by_key('1SacBnsUW4fzqGYl0k5FVV5AFq2TJvUlBOIWp1NzLh88')
book= gc.open('maize_studies')
sheet = book.worksheet("all")
df = pd.DataFrame(sheet.get_all_records())
print(df)

