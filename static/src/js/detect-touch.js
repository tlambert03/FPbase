window.USER_IS_TOUCHING = false

window.addEventListener('touchstart', function onFirstTouch() {
  // we could use a class
  document.body.classList.add('user-is-touching');
  // or set some global variable
  window.USER_IS_TOUCHING = true;
  // we only need to know once that a human touched the screen, so we can stop listening now
  window.removeEventListener('touchstart', onFirstTouch, false);
}, false);
