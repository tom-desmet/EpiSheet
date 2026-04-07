(function () {
  var tip = null;

  function ensureTip() {
    if (tip) return;
    tip = document.createElement('div');
    tip.id = 'epi-tip';
    document.body.appendChild(tip);
  }

  document.addEventListener('click', function (e) {
    ensureTip();
    if (e.target.classList.contains('info-mark')) {
      var r = e.target.getBoundingClientRect();
      tip.textContent = e.target.getAttribute('data-tip');
      tip.style.display = 'block';
      var left = r.left + window.scrollX;
      if (left + 200 > window.innerWidth) left = window.innerWidth - 208;
      tip.style.left = left + 'px';
      tip.style.top  = (r.bottom + window.scrollY + 6) + 'px';
      e.stopPropagation();
    } else {
      tip.style.display = 'none';
    }
  });
})();
