/* Remove the leading blank line in code blocks caused by <pre>\n<code> HTML pattern */
document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('pre > code').forEach(function (code) {
    var prev = code.previousSibling;
    if (prev && prev.nodeType === Node.TEXT_NODE) {
      prev.remove();
    }
  });
});
