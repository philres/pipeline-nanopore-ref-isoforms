#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import sys
import requests as rq

# Parse command line arguments:
parser = argparse.ArgumentParser(description="""Generate an HMM logo using the SkyLign REST API.""")
parser.add_argument('-r', metavar='report_png', type=str, help="Report PDF (skylign_logo.png).", default="SkyLign_logo.png")
parser.add_argument('input', metavar='input_hmm', type=str, help="Input profile HMM.")


if __name__ == '__main__':
    args = parser.parse_args()
    URL = "http://skylign.org"
    cont = open(args.input, 'rb').read()
    files = {'file': cont}
    data = [("processing", "hmm")]
    head = {"Accept": "application/json"}
    print("Submiting profile HMM to SkyLign:", args.input, file=sys.stderr)
    res = rq.post(URL, files=files, data=data, headers=head)
    if res.status_code != rq.codes.ok:
        print("Failed to submit HMM to SkyLign: ", res.content, file=sys.stderr)
        sys.exit(1)
    png_url = json.loads(res.content)["url"]
    png_res = rq.get(png_url, headers={"Accept": "image/png"})
    if png_res.status_code != rq.codes.ok:
        print("Failed to download HMM logo to SkyLign: ", png_res.content, file=sys.stderr)
        sys.exit(1)
    with open(args.r, 'wb') as fh:
        fh.write(png_res.content)
    print("HMM logo written to:", args.r, file=sys.stderr)
