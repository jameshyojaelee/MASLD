"""Explicit descriptive annotation sources; prediction inputs are unaffected."""
import json


def annotation_paths(candidate):
    path = candidate / 'annotation_sources.json'
    config = json.loads(path.read_text()) if path.exists() else {}
    return {key: candidate / config.get(key, key) for key in ['disease_links', 'motif_context']}
