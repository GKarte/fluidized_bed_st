# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 18:51:55 2023

@author: Gregor Karte
"""

from io import BytesIO
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd

def fig2img(fig):
    """
    Konvertiert matplotlib figure in ein png.

    Parameters
    ----------
    fig : matplotlib.figure
        Matplotlib Plot.

    Returns
    -------
    img : png
        png des plots.

    """
    buf = BytesIO()
    fig.savefig(buf, format='png')
    img = buf
    img.seek(0)
    # img = Image.open(buf)
    return img

def str_date_time():
    # using now() to get current time
    current_time = datetime.datetime.now()
    string = f'{current_time.day}{current_time.month}{current_time.year}_{current_time.hour}{current_time.minute}'
    return string


def read_ipse_txt(txt_file):
    df = pd.read_csv(txt_file, sep="\t\t", header=None, index_col=0, engine='python')
    return df


if __name__ == "__main__":
    
    print(read_ipse_txt(r"C:\Users\Gregor\Downloads\PG_SER_20MW.txt"))