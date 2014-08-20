# timer.py
# -------------------------
# May 31, 2013; Alex Safatli
# -------------------------
# Time Estimation Library

import time

def getTime():
    return time.time()

def estimateRemainingTime(starttime,numAt,totalNum):
    if numAt != 0: return ((time.time() - starttime)/(numAt))*(totalNum-numAt)
    else: return 0