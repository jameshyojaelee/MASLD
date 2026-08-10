/**
 * 对目标className的元素增加CDA气泡呼出事件。
 * 气泡在点击其他非目标className的元素或者滚动时隐藏。
 * @param {string} targetClassName 目标能呼出此气泡的、特有的className。
 */
// eslint-disable-next-line
function addHMSMAPopup(targetClassName, type) {
  var className = targetClassName || 'file-name';

  // 创建隐藏的气泡元素。
  var popup = document.createElement('div');
  popup.id = 'FileNamePopup';
  popup.className = 'q-menu q-position-engine';
  var popupCard = document.createElement('div');
  popupCard.className = 'q-card';
  var popupCardSection = document.createElement('div');
  popupCardSection.className = 'q-card__section q-card__section--vert';
  popup.appendChild(popupCard);
  popupCard.appendChild(popupCardSection);
  if (type) {
    popupCardSection.innerHTML = 'The data is requested to be embargoed until June 30, 2026.';
  } else
    popupCardSection.innerHTML =
      'HMSMA database is under control. Please contact <a href="mailto:jin.chai@cldcsw.org" target="_blank">jin.chai@cldcsw.org</a>.';
  document.body.appendChild(popup);

  // 点击目标className元素时呼出气泡，否则隐藏气泡。
  window.addEventListener('click', function (event) {
    if (event.target.className.indexOf(className) >= 0) {
      popup.style = 'visibility: visible; top: ' + event.clientY + 'px; left: ' + event.clientX + 'px;';
      if (className === 'data-name') {
        popup.style = 'visibility: visible; top: ' + event.clientY + 'px; left: ' + (event.clientX - 380) + 'px;';
      }
    } else popup.style = '';
  });

  // 滚动时隐藏气泡。
  window.addEventListener('scroll', function () {
    popup.style = '';
  });
}
