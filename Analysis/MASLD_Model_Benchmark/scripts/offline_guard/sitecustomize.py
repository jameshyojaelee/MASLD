"""Disable IPv4/IPv6 sockets for offline benchmark model execution."""

from __future__ import annotations

import socket


_SOCKET = socket.socket


class GuardedSocket(_SOCKET):
    def __new__(cls, family: int = -1, *args: object, **kwargs: object) -> object:
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for benchmark model execution")
        return super().__new__(cls, family, *args, **kwargs)


socket.socket = GuardedSocket


def blocked_connection(*args: object, **kwargs: object) -> object:
    raise RuntimeError("network disabled for benchmark model execution")


socket.create_connection = blocked_connection
