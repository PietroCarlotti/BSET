(function () {
  function applyFixes() {

    /* Fix 1 — Remove blank first line in code blocks.
       downlit generates <pre>\n<code>, and the \n renders as a visible
       empty line before the first line of code. */
    document.querySelectorAll('pre > code').forEach(function (code) {
      var prev = code.previousSibling;
      if (prev && prev.nodeType === Node.TEXT_NODE) {
        prev.remove();
      }
    });

    /* Fix 2 — Center tables.
       Bootstrap sets width:100% on .table which makes margin:auto ineffective.
       Wrapping in a flex container is the only reliable fix. */
    document.querySelectorAll('table.table').forEach(function (table) {
      if (!table.parentNode || table.parentNode.dataset.tableWrapper) return;
      table.style.width = 'auto';  // override Bootstrap's width:100%
      var wrapper = document.createElement('div');
      wrapper.dataset.tableWrapper = '1';
      wrapper.style.cssText = 'display:flex;justify-content:center;width:100%;';
      table.parentNode.insertBefore(wrapper, table);
      wrapper.appendChild(table);
    });

  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', applyFixes);
  } else {
    applyFixes();
  }
})();
