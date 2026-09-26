"""Command line entry points; no CFX, GPU, or commercial software needed."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .workflows import reproduce, sweep, validate_cfd, paired_comparison


def main(argv=None):
    parser=argparse.ArgumentParser(description="K1 public-data reproduction and purge-allocation screening")
    commands=parser.add_subparsers(dest="command",required=True)
    for name in ("reproduce","sweep","demo","validate-cfd","paired"):
        p=commands.add_parser(name)
        p.add_argument("--out",type=Path,required=True,help="New/empty output directory")
        if name in ("reproduce","sweep","demo"):
            p.add_argument("--no-plots",action="store_true")
        if name in ("sweep","demo"):
            p.add_argument("--baseline",type=Path,required=True,help="Explicit baseline JSON; illustrative values stay labeled")
            p.add_argument("--budget-scale",type=float,default=1.0)
            p.add_argument("--points",type=int,default=201)
        if name in ("validate-cfd","paired"):
            p.add_argument("--results",type=Path,required=True)
        if name=="paired":
            p.add_argument("--reference-allocation",required=True)
    args=parser.parse_args(argv)
    try:
        if args.command=="reproduce":
            result=reproduce(args.out,not args.no_plots)
        elif args.command=="sweep":
            result=sweep(args.baseline,args.out,args.budget_scale,args.points,not args.no_plots)
        elif args.command=="demo":
            if args.out.exists() and any(args.out.iterdir()):
                raise ValueError("Demo output directory must be new or empty.")
            result={"fit":reproduce(args.out/"reproduction",not args.no_plots),
                    "sweep":sweep(args.baseline,args.out/"allocation",args.budget_scale,args.points,not args.no_plots)}
        elif args.command=="validate-cfd":
            result=validate_cfd(args.results,args.out)
        else:
            result=paired_comparison(args.results,args.reference_allocation,args.out)
    except (ValueError,KeyError,OSError,TypeError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        return 2
    print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
