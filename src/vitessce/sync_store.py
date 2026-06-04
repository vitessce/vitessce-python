# Copied from zarr-python
# Reference: https://github.com/zarr-developers/zarr-python/blob/fe229107f9915f05817f7a664d3550695ff9ca44/src/zarr/testing/stateful.py#L438

import builtins
from typing import Any
import zarr
from zarr.abc.store import Store
from zarr.core.buffer import Buffer, BufferPrototype


class SyncStoreWrapper(zarr.core.sync.SyncMixin):
    def __init__(self, store: Store) -> None:
        """Synchronous Store wrapper

        This class holds synchronous methods that map to async methods of Store classes.
        The synchronous wrapper is needed because hypothesis' stateful testing infra does
        not support asyncio so we redefine sync versions of the Store API.
        https://github.com/HypothesisWorks/hypothesis/issues/3712#issuecomment-1668999041
        """
        self.store = store

    @property
    def read_only(self) -> bool:
        return self.store.read_only

    def set(self, key: str, data_buffer: Buffer) -> None:
        return self._sync(self.store.set(key, data_buffer))

    def list(self) -> builtins.list[str]:
        return self._sync_iter(self.store.list())

    def get(self, key: str, prototype: BufferPrototype, **kwargs) -> Buffer | None:
        return self._sync(self.store.get(key, prototype=prototype, **kwargs))

    def get_partial_values(
        self, key_ranges: builtins.list[Any], prototype: BufferPrototype
    ) -> builtins.list[Buffer | None]:
        return self._sync(self.store.get_partial_values(prototype=prototype, key_ranges=key_ranges))

    def delete(self, path: str) -> None:
        return self._sync(self.store.delete(path))

    def is_empty(self, prefix: str) -> bool:
        return self._sync(self.store.is_empty(prefix=prefix))

    def clear(self) -> None:
        return self._sync(self.store.clear())

    def exists(self, key: str) -> bool:
        return self._sync(self.store.exists(key))

    def list_dir(self, prefix: str) -> None:
        raise NotImplementedError

    def list_prefix(self, prefix: str) -> None:
        raise NotImplementedError

    @property
    def supports_listing(self) -> bool:
        return self.store.supports_listing

    @property
    def supports_writes(self) -> bool:
        return self.store.supports_writes

    @property
    def supports_deletes(self) -> bool:
        return self.store.supports_deletes
