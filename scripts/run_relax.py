#!/usr/bin/env python3

import argparse
import pyrosetta

from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects


parser = argparse.ArgumentParser()
parser.add_argument("input_pdb")
parser.add_argument("output_pdb")
parser.add_argument("--xml", required=True)
args = parser.parse_args()

pyrosetta.init("-corrections::beta_nov16 true")

pose = pyrosetta.pose_from_pdb(args.input_pdb)

xml_objects = XmlObjects.create_from_file(args.xml)
protocol = xml_objects.get_mover("ParsedProtocol")

protocol.apply(pose)

pose.dump_pdb(args.output_pdb)

print(f"Wrote {args.output_pdb}")